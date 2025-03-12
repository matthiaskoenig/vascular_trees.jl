module Simulation_Helpers

export create_simulations, create_benchmarked_simulations

using OrdinaryDiffEq, DataFrames, CSV, Sundials
using TimerOutputs, Plots

include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR, MODEL_PATH 
import .Utils.Definitions: flow_directions, tree_definitions, ODE_groups
import .Utils.Benchmarking: save_times_as_csv

# using ..Utils: JULIA_RESULTS_DIR, MODEL_PATH
# using ..Utils.Definitions: flow_directions, ODE_groups
# using ..Utils.Benchmarking: save_times_as_csv

include("../models/julia_from_jgraph.jl")
import .Julia_from_jgraph: get_ODE_components

include("../../" * MODEL_PATH)
import .Julia_models: jf_dxdt!

# Already specified in utils.jl
const flow_direction = flow_directions()
const trees = tree_definitions()
const groups = ODE_groups()
const tspan = Vector{Float64}(undef, 2)

function solve_ODE(x0, p, tspan, sol_options)

    problem = ODEProblem(jf_dxdt!, x0, tspan, p)

    solution = solve(
        problem,
        sol_options.solver,
        #saveat = [tspan[2]],
        dense = false,
        reltol = sol_options.relative_tolerance,
        abstol = sol_options.absolute_tolerance,
        #save_idxs = idxs_tosave,
        #progress = true
    )

    #print(solution.destats)
    # print(solution)
    # display(plot(solution))

    return solution

end

function create_simulations(g_options, sim_options, sol_options)
    n_intervals = Integer(sim_options.tspan[2] / sim_options.sdt)
    for n_node ∈ g_options.n_nodes, tree_id ∈ g_options.tree_ids
        graph_id = "$(tree_id)_$(n_node)"
        vascular_trees = values(trees.vascular_trees[tree_id])

        u0 = Dict(vascular_tree => [0.0] for vascular_tree in Iterators.flatten(vascular_trees))
        p = Dict{String, Tuple}()
        synch_idxs = Dict{String, Vector{Integer}}()
        variable_ids = Dict{String, Vector{Symbol}}()
        solutions = Dict(vascular_tree => [DataFrame() for _ in 1:n_intervals+1] for vascular_tree in Iterators.flatten(vascular_trees))
        solutions["T"] = [DataFrame() for _ in 1:n_intervals+1]

        # preparation step - get starting initial values and get parameters for all trees
        for vascular_tree ∈ Iterators.flatten(vascular_trees)
            u0[vascular_tree], p[vascular_tree] = get_ODE_components(tree_id, n_node, vascular_tree)
            synch_idxs[vascular_tree] = get_synchronization_indices(vascular_tree, p[vascular_tree])
            variable_ids[vascular_tree] = get_variable_ids(vascular_tree, p[vascular_tree][3])
        end
        u0_terminal, p_terminal = get_ODE_components(tree_id, n_node, "T")
        variable_ids["T"] = get_variable_ids("T", length(u0_terminal[1, :]))

        solve_trees_separately!(solutions, u0_terminal, p_terminal, u0, p, sim_options, sol_options, vascular_trees, synch_idxs)
        for key in keys(solutions)
            solutions[key][end] = reduce(vcat, solutions[key])
        end
        (size(solutions["A"][end])[2] == size(solutions["T"][end])[2]) ? (print(true)) : (print(false))
        
    end
end

function solve_trees_separately!(solutions, u0_terminal, p_terminal, u0, p, sim_options, sol_options, vascular_trees, synch_idxs)
    tmin = sim_options.tspan[1]
    tmax = sim_options.tspan[2]
    sdt = sim_options.sdt
    ki = 1
    t_left = tmin
    t_right = sdt

    while t_right <= tmax
        tspan .= [t_left, t_right]
        for vascular_tree ∈ Iterators.flatten(vascular_trees)
            solution = solve_ODE(u0[vascular_tree], p[vascular_tree], tspan, sol_options)
            u0[vascular_tree] .= solution[end]
            solutions[vascular_tree][ki] = DataFrame(solution)
        end
        solution_terminal = solve_ODE(u0_terminal, p_terminal, tspan, sol_options)
        u0_terminal .= solution_terminal[end]
        solutions["T"][ki] = DataFrame(solution_terminal)

        for (ki, inflow_tree) in enumerate(p_terminal[2])
            if inflow_tree ∈ flow_direction.inflow_trees
                u0_terminal[ki, :] .= view(u0[inflow_tree], synch_idxs[inflow_tree])
            end
        end

        for outflow_tree in flow_direction.outflow_trees
            u0[outflow_tree][synch_idxs[outflow_tree]] .= u0_terminal[1, :]
        end

        t_left += sdt
        t_right += sdt
        ki += 1
    end
end

function get_variable_ids(vascular_tree_id::String, element_ids::Vector{Tuple{T,T}}) where {T<:Integer}
    variable_ids = [Symbol("$(vascular_tree_id)_$(element_id[1])_$(element_id[2])") for element_id in element_ids]
    return variable_ids
end

function get_variable_ids(vascular_tree_id::String, n_elements::Integer)
    variable_ids = [Symbol("$(vascular_tree_id)_$(kv)") for kv in eachindex(n_elements)]
    return variable_ids
end

function get_synchronization_indices(graph_id, graph_parameters)
    if graph_id in flow_direction.inflow_trees
        synch_idx = get_indices(graph_parameters[6], [groups.preterminal])
    else
        synch_idx = get_indices(graph_parameters[6], [groups.terminal])
    end
    return synch_idx
end

function create_simulations(g_options, sim_options, bench_options, sol_options)
    to = TimerOutput()
    graph_id::String = ""

    for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
        graph_id = "$(tree_id)_$(n_node)"
        vascular_trees = values(trees.vascular_trees[tree_id])

        for vascular_tree ∈ Iterators.flatten(vascular_trees)
            x0, graph_p = get_ODE_components(tree_id, n_node, vascular_tree)
            p = (
                graph_p.is_inflow,
                graph_p.flows,
                graph_p.volumes,
                graph_p.ODE_groups,
                graph_p.pre_elements,
                graph_p.post_elements,
            )
            idxs_tosave = get_indices(p[4], [groups.marginal, groups.terminal])
            @timeit to "$(graph_id)_$(vascular_tree)" begin
                for ki = 1:bench_options.n_iterations
                    @timeit to "$(graph_id)_$ki" ODE_solver(x0, p, sim_options, sol_options, idxs_tosave)
                end
            end
            show(to, sortby = :firstexec)
            if bench_options.save_running_times
                n_species = length(x0)
                n_terminals = Int64(n_species / 3)
                save_times_as_csv(
                    to,
                    n_node,
                    "$(tree_id)_$(vascular_tree)",
                    n_species,
                    sol_options.solver_name,
                    n_terminals,
                )
            end
            reset_timer!(to::TimerOutput)
        end
    end
end

function save_simulations_to_csv(
    simulations,
    column_names::Vector{String},
    simulations_path::String,
)
    # Convert solution to DataFrame
    df = DataFrame(simulations)

    # Write DataFrame to CSV
    header = vcat(["time"], column_names)
    CSV.write(simulations_path, df, header = header)
end

function get_indices(collection, groups)
    indices::Vector{Int64} = []
    @inbounds for (ke, element) ∈ enumerate(collection)
        (element ∈ groups) && (push!(indices, ke))
    end
    return indices
end    
end
