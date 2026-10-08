"""
Functions that set up and run the coupled simulation of the vascular trees and save the results.

Data flow of `run_simulations`:
1. Preparation: ODE parameters, initial values and species ids are loaded for every vascular
   tree (A, P, V, B, ...) and for the terminal part ("T").
2. Integration (`solve_tree!`): time is advanced in steps of `dt`. In every step each tree and
   the terminal part are integrated separately and then coupled through the values at the
   connecting species (see `solve_tree!`).
3. Saving: for each subsystem a CSV file with columns `t` + species ids and one row per time
   point is written (only if `sim_options.save_simulations = true` and no benchmarking).

Subsystem = vascular tree (e.g. "A") or the terminal part ("T").
"""
module Helpers

export run_simulations

using DifferentialEquations, DataFrames, CSV, Dictionaries, TimerOutputs, ProgressMeter

using ..Paths: JULIA_RESULTS_DIR, MODEL_PATH

using ..Definitions: flow_directions, ODE_groups, vascular_tree_parameters
using ..Benchmarking: save_times_as_csv

using ..get_ode_parameters: get_ODE_parameters, get_initial_values

# Which trees are inflow / outflow, and the numbers marking the ODE groups of the species
# (both defined in utils.jl)
const flow_direction = flow_directions()
const groups = ODE_groups()

using ..Pharmacokinetic_models: jf_dxdt!

using Profile
using PProf


"""
    run_simulations(tree_info, sim_options, sol_options, additional_sol_options,
                    flow_scaling_factor, bench_options)

Create and run the coupled tree simulation for one graph and one flow scaling factor.

If `sim_options.benchmark` is `true`, the simulation is run `bench_options.n_iterations` times
and timed (the results are printed and, if `bench_options.save_running_times`, saved to a CSV).
Otherwise it is run once and, if `sim_options.save_simulations`, the solutions are saved to CSV.

# Arguments
- `tree_info`: `Tree_structure` of the graph (ids of the trees, graph directory, ...).
- `sim_options`: simulation options (`tspan`, `steps`, `dt`, `save_simulations`, `benchmark`).
- `sol_options`: ODE solver options (solver, tolerances).
- `additional_sol_options`: additional integrator arguments (not mandatory), a `NamedTuple`.
- `flow_scaling_factor`: factor by which the flow is multiplied, 1.0 = original flow.
- `bench_options`: benchmark options (`n_iterations`, `save_running_times`).
"""
function run_simulations(
    tree_info,
    sim_options,
    sol_options,
    additional_sol_options,
    flow_scaling_factor::AbstractFloat,
    bench_options,
)

    # TODO: mostly all of this can be preallocated

    # graph_subsystem - vascular trees + terminal part ("T")
    graph_subsystems = [tree_info.vascular_trees; ["T"]]
    # solutions[graph_subsystem][k] - values of all species at the k-th time point (first value is t)
    solutions = dictionary(
        graph_subsystem => [Float64[] for _ = 1:sim_options.steps+1] for
        graph_subsystem in graph_subsystems
    )
    species_ids = similar(solutions, Vector{Symbol})

    # preparation step - get parameters and species ids (as symbols) for all trees
    p = Dict{String,vascular_tree_parameters}()
    for vascular_tree ∈ tree_info.vascular_trees
        p[vascular_tree] = get_ODE_parameters(tree_info, vascular_tree, flow_scaling_factor)
        species_ids[vascular_tree] = collect_species_ids(p[vascular_tree].species_ids)
    end
    # starting initial values of trees
    u0 = dictionary(
        vascular_tree => get_initial_values(p[vascular_tree].ODE_groups) for
        vascular_tree in tree_info.vascular_trees
    )
    # indices of the species that connect trees and terminal nodes
    synch_idxs = dictionary(
        vascular_tree => get_synchronization_indices(vascular_tree, p[vascular_tree].ODE_groups)
        for vascular_tree in tree_info.vascular_trees
    )

    # parameters, initial values and species ids for terminal nodes
    # they are stored separately from trees, because they are stored differently
    p_terminal = get_ODE_parameters(tree_info, "T", flow_scaling_factor)
    u0_terminal = get_initial_values(p_terminal.flow_values)
    species_ids["T"] = collect_species_ids(vec(p_terminal.species_ids))

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
        Profile.Allocs.clear()
        Profile.Allocs.@profile sample_rate = 0.0001 solve_tree!(
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
        PProf.Allocs.pprof(; from_c=false)
        if sim_options.save_simulations
            save_solutions(solutions, species_ids, tree_info, sim_options.dt, flow_scaling_factor)
        end
    end
end


"""
    solve_tree!(solutions, u0_terminal, p_terminal, u0, p, sim_options, sol_options,
                synch_idxs, additional_sol_options)

Integrate the coupled system of vascular trees and terminal part over `sim_options.tspan`.

Every time step of length `dt`:
1. each tree and the terminal part are integrated separately from their current values;
2. the values at the start of the step are stored in `solutions` (with `t` as first value);
3. the coupling values are exchanged for the next step: the preterminal species of the
   inflow trees (A, P) become the inflow of the terminal part, and the first row
   of the terminal part (its outflow) becomes the values of the terminal species of
   the outflow trees (V, B).

`solutions`, `u0` and `u0_terminal` are modified in place.

# Arguments
- `solutions`: for every subsystem a vector (length `steps+1`) to store the values of each time point.
- `u0_terminal`, `p_terminal`: initial values and parameters of the terminal part.
- `u0`, `p`: initial values and parameters of the trees, with tree ids as keys.
- `sim_options`, `sol_options`, `additional_sol_options`: see [`run_simulations`](@ref).
- `synch_idxs`: indices of the species that connect each tree with the terminal part
  (see [`get_synchronization_indices`](@ref)).
"""
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

    progress = Progress(
        Int(sim_options.steps);
        dt = dt,
        barglyphs = BarGlyphs('|', '█', ['▁', '▂', '▃', '▄', '▅', '▆', '▇'], ' ', '|'),
        color = :magenta,
    )
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
        # update in terminal part (inflow from the inflow trees)
        for (ki, species_id) in enumerate(view(p_terminal.species_ids, :, 1))
            vascular_tree_id = first(species_id, 1)
            if vascular_tree_id ∈ flow_direction.inflow_trees
                u0_terminal[ki, :] .=
                    view(u0[vascular_tree_id], synch_idxs[vascular_tree_id])
            end
        end
        # update for outflow trees (outflow from the terminal part)
        for outflow_tree in flow_direction.outflow_trees
            u0[outflow_tree][synch_idxs[outflow_tree]] .= view(u0_terminal, 1, :)
        end

        # updating state for next integration step
        t = t + dt
        kl += 1
        next!(progress)
    end
end


"""
    collect_species_ids(species_ids)

Convert species ids (strings) to symbols and add `:t` (time) as the first one, so that the
result can be used as column names of the simulation table.
"""
function collect_species_ids(species_ids::Array{String})
    return [:t; Symbol.(species_ids)]
end


"""
    get_synchronization_indices(graph_id, species_ODE_groups)

Indices of the species of the tree `graph_id` that connect it with the terminal part:
preterminal species for inflow trees, terminal species for outflow trees.
`species_ODE_groups` holds the ODE group of every species of the tree.
"""
function get_synchronization_indices(graph_id, species_ODE_groups)
    if graph_id in flow_direction.inflow_trees
        return get_indices(species_ODE_groups, [groups.preterminal])
    else
        return get_indices(species_ODE_groups, [groups.terminal])
    end
end


"""
    get_indices(collection, targets)

Indices of the elements of `collection` that are in `targets`.
"""
function get_indices(collection, targets)
    return findall(in(targets), collection)
end


"""
    simulation_csv_path(tree_info, graph_subsystem, flow_scaling_factor, dt)

Path of the CSV file with the simulation of `graph_subsystem`. The flow scaling factor
is part of the file name only if it differs from 1.0.
"""
function simulation_csv_path(tree_info, graph_subsystem, flow_scaling_factor, dt)
    flow_part = flow_scaling_factor == 1.0 ? "" : "_$(flow_scaling_factor)Q"
    return joinpath(
        tree_info.GRAPH_DIR,
        "simulations",
        "$(graph_subsystem)$(flow_part)_simulations_$(dt)_dt.csv",
    )
end


"""
    save_solutions(solutions, species_ids, tree_info, dt, flow_scaling_factor)

Save the solution of every subsystem as a CSV table (columns: `species_ids`, one row per
time point), see [`simulation_csv_path`](@ref) for the file names.
"""
function save_solutions(solutions, species_ids, tree_info, dt, flow_scaling_factor)
    @info "Saving results"
    for (graph_subsystem, solution) in pairs(solutions)
        # skip time points that were not filled
        filter!(!isempty, solution)
        # turn the vector of states into columns: one column per species
        columns = Vector{Vector}(undef, 0)
        for species_idx in eachindex(first(solution))
            push!(columns, getindex.(solution, species_idx))
        end
        table = DataFrame(columns, species_ids[graph_subsystem], copycols = false)
        save_simulations_to_csv(
            table,
            simulation_csv_path(tree_info, graph_subsystem, flow_scaling_factor, dt),
        )
    end
end


"""
    save_simulations_to_csv(simulations, simulations_path)

Write the `simulations` table to the CSV file `simulations_path`.
"""
function save_simulations_to_csv(simulations::DataFrame, simulations_path::String)
    CSV.write(simulations_path, simulations)
end

end
