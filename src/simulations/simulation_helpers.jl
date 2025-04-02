module Simulation_Helpers

export run_simulations, create_benchmarked_simulations

using DifferentialEquations,
    DataFrames, CSV, Sundials, Dictionaries, TimerOutputs

using ..Utils: MODEL_PATH
using ..Utils.Definitions: flow_directions, ODE_groups, terminal_parameters, vascular_tree_parameters
using ..Utils.Benchmarking: save_times_as_csv

using ..Julia_from_jgraph: get_ODE_parameters, get_initial_values

# Already specified in utils.jl
const flow_direction = flow_directions()
const groups = ODE_groups()
const tspan = Vector{Float64}(undef, 2)

include("../../" * MODEL_PATH)
using .Pharmacokinetic_models: jf_dxdt!

using InteractiveUtils

function run_simulations(
    tree_info,
    sim_options,
    sol_options,
    additional_sol_options,
    flow_scaling_factor::AbstractFloat,
    bench_options,
)
    """Create and run coupled tree simulations.

    t_options: Graph options
    sim_options: Simulation options
    sol_options: ODE Solver options
    additional_sol_options:: Additional integrator arguments which are not mandatory
    bench_options:: Benchmark options
    """

    # MOSTLY ALL OF THIS STUFF CAN BE PREALLOCATED
    
    #u0 = dictionary(vascular_tree => [0.0] for vascular_tree in tree_info.vascular_trees)
    # create dictionary to store the parameters of trees (placeholder)
    p = Dict{String, vascular_tree_parameters}() #dictionary(vascular_tree => vascular_tree_parameters() for vascular_tree in tree_info.vascular_trees)

    # create dictionaries to store solutions and species ids
    # graph_subsystem - vascular trees + terminal part
    graph_subsystems = [tree_info.vascular_trees; ["T"]]
    solutions = dictionary(
        graph_subsystem => [Float64[] for _ = 1:sim_options.steps+1] for
        graph_subsystem in graph_subsystems
    )
    species_ids = similar(solutions, Vector{Symbol})

    # preparation step - get starting initial values, 
    # parameters and species ids (as symbols) for all trees
    for vascular_tree ∈ tree_info.vascular_trees
        p[vascular_tree] = get_ODE_parameters(tree_info, vascular_tree, flow_scaling_factor)
        species_ids[vascular_tree] = collect_species_ids(p[vascular_tree].species_ids) # p[vascular_tree][3] - Vector with String species ids
    end
    # create dictionary that stores the initial values of trees (placeholder)
    u0 = dictionary(vascular_tree => get_initial_values(p[vascular_tree].ODE_groups) for vascular_tree in tree_info.vascular_trees) # p[vascular_tree][3] - Vector with String species ids
    # create dictionary that stores the indices of the species 
    # that connect trees and terminal nodes
    synch_idxs = dictionary(vascular_tree => get_synchronization_indices(vascular_tree, p[vascular_tree].ODE_groups) for vascular_tree in tree_info.vascular_trees)

    # get starting initial values, parameters and
    # species ids (as symbols) for terminal nodes
    # u0 and p for terminal nodes are stored as separate values (not with trees),
    # because they are stored differently 
    p_terminal = get_ODE_parameters(tree_info, "T", flow_scaling_factor)
    u0_terminal = get_initial_values(p_terminal.flow_values)

    species_ids["T"] = collect_species_ids(vec(p_terminal.x_affiliations)) # p_terminal[2] - Matrix with String species ids

    if sim_options.benchmark
        # benchmark solving function
        to::TimerOutput = TimerOutput()
        @timeit to tree_info.graph_id begin
            for ki = 1:bench_options.n_iterations
                @timeit to "$(tree_info.graph_id)_$ki" solve_tree!(
                    solutions,
                    u0_terminal,
                    p_terminal,
                    u0,
                    p,
                    sim_options,
                    sol_options,
                    synch_idxs,
                    additional_sol_options,
                )
            end
        end

        show(to, sortby = :firstexec)
        if bench_options.save_running_times
            @info "Pay attention to the number of species, they can be calculated wrong"
            n_terminal::Integer = size(u0_terminal)[1]
            save_times_as_csv(
                to,
                tree_info.n_node,
                tree_info.graph_id,
                n_terminal * 9,
                sol_options.solver_name,
                n_terminal,
            )
        end
        reset_timer!(to)
    else
        # solve the tree
        solve_tree!(
            solutions,
            u0_terminal,
            p_terminal,
            u0,
            p,
            sim_options,
            sol_options,
            synch_idxs,
            additional_sol_options,
        )
        if sim_options.save_simulations
            @info "Saving results"
            for (graph_subsystem, solution) in pairs(solutions)
                filter!(!isempty, solution)
                columns = Vector{Vector}(undef, 0)
                for species_idx in eachindex(first(solution))
                    push!(columns, getindex.(solution, species_idx))
                end
                solution = DataFrame(columns, species_ids[graph_subsystem], copycols=false)
                #rename!(solution, species_ids[graph_subsystem])
                if flow_scaling_factor == 1.0
                    path = joinpath(
                        tree_info.GRAPH_DIR,
                        "simulations",
                        "$(graph_subsystem)_simulations_$(sim_options.dt)_dt.csv",
                    )
                else
                    path = joinpath(
                        tree_info.GRAPH_DIR,
                        "simulations",
                        "$(graph_subsystem)_$(flow_scaling_factor)Q_simulations_$(sim_options.dt)_dt.csv",
                    )
                end
                save_simulations_to_csv(
                    solution,
                    path
                )
            end
        end
    end
end

function solve_tree!(
    solutions,
    u0_terminal,
    p_terminal,
    u0,
    p,
    sim_options,
    sol_options,
    synch_idxs,
    additional_sol_options,
)

    # storing integrators for reuse
    integrators = Dict{String,Any}()
    # collecting integrator arguments together
    integrator_options::NamedTuple = (
        reltol = sol_options.relative_tolerance,
        abstol = sol_options.absolute_tolerance,
        save_end = false,
        additional_sol_options...
    )

    # setup integrators for tree problems
    for vascular_tree_id ∈ keys(u0)
        problem = ODEProblem(jf_dxdt!, u0[vascular_tree_id], [0, 0.01], p[vascular_tree_id])
        integrators[vascular_tree_id] =
            init(problem, sol_options.solver; integrator_options...)
    end

    # setup integrator for terminal node problem
    problem_terminal = ODEProblem(jf_dxdt!, u0_terminal, [0, 0.01], p_terminal)
    integrator_terminal = init(problem_terminal, sol_options.solver; integrator_options...)

    @info "Running integration"
    tmin = sim_options.tspan[1]
    tmax = sim_options.tspan[2]
    dt = sim_options.dt
    kl = 1  # loop iterator
    t = tmin

    while t <= tmax

        # solve current step
        for (vascular_tree_id, integrator) ∈ pairs(integrators)
            # reinitialize trees' integrator problem
            reinit!(
                integrator,
                u0[vascular_tree_id];
                t0 = t,
                tf = t + dt,
                erase_sol = true
            )
            # solve the timestep for the subtree
            solve!(integrator)
            # final values at end of integration
            u0[vascular_tree_id] .= integrator.uprev
        end
        # reinitialize terminal integrator problem
        reinit!(integrator_terminal, u0_terminal; t0 = t, tf = t + dt, erase_sol = true)
        # solve the timestep for the terminal nodes
        solve!(integrator_terminal)
        # final values at end of integration
        u0_terminal .= integrator_terminal.uprev

        # store numerical solutions
        for (vascular_tree_id, integrator) ∈ pairs(integrators)
            solutions[vascular_tree_id][kl] = [t ; Vector(integrator.sol)]
        end
        solutions["T"][kl] = [t ; vec(integrator_terminal.sol)]

        # updating initial values
        # update in terminal part
        for (ki, species_id) in enumerate(view(p_terminal.x_affiliations, :, 1))
            vascular_tree_id = first(species_id, 1)
            if vascular_tree_id ∈ flow_direction.inflow_trees
                u0_terminal[ki, :] .=
                    view(u0[vascular_tree_id], synch_idxs[vascular_tree_id])
            end
        end
        # update for outflow trees
        for outflow_tree in flow_direction.outflow_trees
            u0[outflow_tree][synch_idxs[outflow_tree]] .= view(u0_terminal, 1, :)
        end

        # updating state for next integration step
        t = t + dt
        kl += 1
    end
end

function collect_species_ids(species_ids::Array{String})
    species_ids = Symbol.(species_ids)
    #species_ids = [Symbol(element_id) for element_id in species_ids]
    pushfirst!(species_ids, :t)
    return species_ids
end

function get_synchronization_indices(graph_id, species_ODE_groups)
    if graph_id in flow_direction.inflow_trees
        synch_idx = get_indices(species_ODE_groups, [groups.preterminal])
    else
        synch_idx = get_indices(species_ODE_groups, [groups.terminal])
    end
    return synch_idx
end

function get_indices(collection, groups)
    indices::Vector{Int64} = []
    @inbounds for (ke, element) ∈ enumerate(collection)
        (element ∈ groups) && (push!(indices, ke))
    end
    return indices
end

function save_simulations_to_csv(simulations::DataFrame, simulations_path::String)
    # Write DataFrame to CSV
    CSV.write(simulations_path, simulations)
end

end
