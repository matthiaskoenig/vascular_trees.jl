"""
Functions that set up and run the coupled simulation of the vascular trees and save the results.

Data flow of `run_simulations`:
1. Preparation: ODE parameters, initial values and species ids are loaded for every vascular
   tree (A, P, V, B, ...) and for the terminal part ("T").
2. Integration (`solve_tree!`): time is advanced in steps of `dt`. In every step each tree and
   the terminal part are integrated separately and then coupled through the values at the
   connecting species (see `solve_tree!`).
3. Saving: for each subsystem an ARROW file with one row per species (column `species_ids`) and
   one column per time point (`t_<time>`) is written (only if `sim_options.save_simulations = true`
   and no benchmarking).

Subsystem = vascular tree (e.g. "A") or the terminal part ("T").
"""
module SimulationHelpers

    export run_simulations

    using DifferentialEquations, DataFrames, CSV, Dictionaries, TimerOutputs, ProgressMeter, Arrow

    using ..Paths: JULIA_RESULTS_DIR, MODEL_PATH

    using ..Definitions: FlowDirections, ODE_groups, VascularTreeParameters
    using ..Benchmarking: save_times_as_csv
    using ..TableUtils: save_as_arrow

    using ..ODEparameters: get_ODE_parameters, get_initial_values

    # Which trees are inflow / outflow, and the numbers marking the ODE groups of the species
    # (both defined in utils.jl)
    const flow_direction = FlowDirections()
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
    - `tree_info`: `TreeStructure` of the graph (ids of the trees, graph directory, ...).
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
        # 1. parameters of trees and of the terminal part
        p = dictionary(
            vascular_tree => get_ODE_parameters(tree_info, vascular_tree, flow_scaling_factor) for
            vascular_tree in tree_info.vascular_trees
        )
        # terminal nodes are stored separately from trees, because they are stored differently
        p_terminal = get_ODE_parameters(tree_info, "T", flow_scaling_factor)

        # 2. initial values and synchronization indices
        u0 = dictionary(
            vascular_tree => get_initial_values(p[vascular_tree].ODE_groups) for
            vascular_tree in tree_info.vascular_trees
        )
        u0_terminal = get_initial_values(p_terminal.flow_values)
        # indices of the species that connect trees and terminal nodes
        synch_idxs = dictionary(
            vascular_tree => get_synchronization_indices(vascular_tree, p[vascular_tree].ODE_groups)
            for vascular_tree in tree_info.vascular_trees
        )

        # 3. preallocated solutions
        # solutions[graph_subsystem][k] - values of all species at the k-th time point
        # (one vector per time point: recording a step is one contiguous copy, and Arrow writes
        # plain Vectors without copying them)
        n_time_points::Int = sim_options.steps + 1
        solutions = dictionary(
            tree => [zeros(length(u0[tree])) for _ = 1:n_time_points] for tree in tree_info.vascular_trees
        )
        insert!(solutions, "T", [zeros(length(u0_terminal)) for _ = 1:n_time_points])

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
                @info "Pay attention to the number of species, they could be calculated in a wrong way"
                n_terminal::Int = size(u0_terminal)[1]
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
                # 4. species ids — same keys as solutions; terminal ids in linear (column-major)
                #    order, the order in which the terminal Matrix u is stored in solutions
                species_ids = similar(solutions, Vector{String})
                for vascular_tree ∈ tree_info.vascular_trees
                    species_ids[vascular_tree] = p[vascular_tree].species_ids
                end
                species_ids["T"] = vec(p_terminal.species_ids)
                times = sim_options.tspan[1] .+ (0:n_time_points-1) .* sim_options.dt
                save_solutions(solutions, species_ids, times, tree_info, sim_options.dt, flow_scaling_factor)
            end
        end
    end


    """
        solve_tree!(solutions, u0_terminal, p_terminal, u0, p, sim_options, sol_options,
                    synch_idxs, additional_sol_options)

    Integrate the coupled system of vascular trees and terminal part over `sim_options.tspan`.

    Every time step of length `dt`:
    1. each tree and the terminal part are integrated separately from their current values;
    2. the values at the start of the step are stored in `solutions`;
    3. the coupling values are exchanged for the next step: the preterminal species of the
    inflow trees (A, P) become the inflow of the terminal part, and the first row
    of the terminal part (its outflow) becomes the values of the terminal species of
    the outflow trees (V, B).

    `solutions`, `u0` and `u0_terminal` are modified in place.

    # Arguments
    - `solutions`: for every subsystem `steps+1` vectors, one per time point, each holding the
    values of all species of that subsystem.
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

        # collecting integrator arguments together
        integrator_options = (
            reltol = sol_options.relative_tolerance,
            abstol = sol_options.absolute_tolerance,
            save_start = false,
            save_end = false,
            maxiters = (sim_options.steps + 1) * 100_000,
            additional_sol_options...
        )
        
        @info "Running integration"
        tmin = sim_options.tspan[1]
        tmax = sim_options.tspan[2]
        dt = sim_options.dt
        # integrate over the whole time span; the loop advances it in steps of dt
        # (tmax + 2dt as margin so that float drift in t never steps past tf)
        tspan = (tmin, tmax + 2dt)

        # setup integrators for tree problems
        integrators = dictionary(
            vascular_tree_id => init(
                ODEProblem(jf_dxdt!, u0[vascular_tree_id], tspan, p[vascular_tree_id]),
                sol_options.solver;
                integrator_options...,
            ) for vascular_tree_id in keys(u0)
        )
        # setup integrator for terminal node problem
        problem_terminal = ODEProblem(jf_dxdt!, u0_terminal, tspan, p_terminal)
        integrator_terminal = init(problem_terminal, sol_options.solver; integrator_options...)

        ids = collect(keys(integrators))
        tree_bindings = [(integ = integrators[id], sol = solutions[id], sidx = synch_idxs[id],
                        is_outflow = id in flow_direction.outflow_trees) for id in ids]
        # terminal rows that receive inflow from an inflow tree: (row index, tree id)
        # computed once, because the assignment of terminal rows to trees does not change
        inflow_rows = [(ki, findfirst(==(tid), ids)) for (ki, tid) in
               enumerate(first.(view(p_terminal.species_ids, :, 1), 1)) if
               tid in flow_direction.inflow_trees]   # binding index instead of String

        progress = Progress(
            Int(sim_options.steps);
            dt = 1.0,
            barglyphs = BarGlyphs('|', '█', ['▁', '▂', '▃', '▄', '▅', '▆', '▇'], ' ', '|'),
            color = :magenta,
            enabled = false
        )
        integrate!(solutions["T"], integrator_terminal, inflow_rows, tree_bindings, dt, progress)

    end

    """
        integrate!(solution_terminal, integrator_terminal, inflow_rows, tree_bindings, dt, progress)

    Time loop of [`solve_tree!`](@ref): store the state of every subsystem at each time point
    (`solution_terminal` and `b.sol` of each tree binding, one vector per time point), then
    advance all subsystems by `dt` and exchange the coupling values. The number of time points
    is `length(solution_terminal)`; after the last one is stored, no further step is integrated.
    """
    function integrate!(solution_terminal, integrator_terminal, inflow_rows, tree_bindings, dt, progress)
        n_time_points = length(solution_terminal)
        for kl in 1:n_time_points
            # store numerical solutions (state at the kl-th time point)
            for b in tree_bindings
                copyto!(b.sol[kl], b.integ.u)
            end
            # linear (column-major) order of the terminal Matrix u
            copyto!(solution_terminal[kl], integrator_terminal.u)
            # last time point stored: nothing left to integrate
            kl == n_time_points && break

            # solve current step: advance every subsystem exactly to t + dt
            for b in tree_bindings
                step!(b.integ, dt, true)
            end
            step!(integrator_terminal, dt, true)

            # updating values at the connecting species
            # update in terminal part (inflow from the inflow trees)
            for (ki, bi) in inflow_rows
                b = tree_bindings[bi]
                copyto!(view(integrator_terminal.u, ki, :), view(b.integ.u, b.sidx))
            end
            # update for outflow trees (outflow from the terminal part)
            for b in tree_bindings
                if b.is_outflow
                    b.integ.u[b.sidx] .= view(integrator_terminal.u, 1, :)
                end
            end
            # u was changed from outside the solver: reset its cached derivative
            for b in tree_bindings
                derivative_discontinuity!(b.integ, true)
            end
            derivative_discontinuity!(integrator_terminal, true)

            next!(progress)
        end
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
        simulation_path(tree_info, graph_subsystem, flow_scaling_factor, dt, format)

    Path of the ARROW file with the simulation of `graph_subsystem`. The flow scaling factor
    is part of the file name only if it differs from 1.0.
    """
    function simulation_path(tree_info, graph_subsystem, flow_scaling_factor, dt, format::String)
        flow_part = flow_scaling_factor == 1.0 ? "" : "_$(flow_scaling_factor)Q"
        return joinpath(
            tree_info.GRAPH_DIR,
            "simulations",
            "$(graph_subsystem)$(flow_part)_simulations_$(dt)_dt.$(format)",
        )
    end


    """
        save_solutions(solutions, species_ids, times, tree_info, dt, flow_scaling_factor)

    Save the solution of every subsystem as an ARROW table with one row per species: column
    `species_ids`, then one column per time point named `t_<time>` (`times[k]` for
    `solutions[graph_subsystem][k]`). See [`simulation_path`](@ref) for the file names.

    The time-point vectors are passed to Arrow as they are (`copycols = false`); plain `Vector`s
    are written without an intermediate copy.
    """
    function save_solutions(solutions, species_ids, times, tree_info, dt, flow_scaling_factor)
        @info "Saving results"
        # rounded, so that float noise of tmin + k*dt does not end up in the column names
        column_names = ["species_ids"; ["t_$(round(t; digits = 10))" for t in times]]
        for (graph_subsystem, solution) in pairs(solutions)
            table = DataFrame(
                AbstractVector[species_ids[graph_subsystem], solution...],
                column_names;
                copycols = false,
            )
            save_as_arrow(
                table,
                simulation_path(tree_info, graph_subsystem, flow_scaling_factor, dt, "arrow"),
            )
            # save_simulations_to_csv(
            #     table,
            #     simulation_csv_path(tree_info, graph_subsystem, flow_scaling_factor, dt, "csv"),
            # )
        end
    end


    """
        save_simulations_to_csv(simulations, simulations_path)

    Write the `simulations` table to the CSV file `simulations_path`.
    """
    function save_simulations_to_csv(simulations::DataFrame, simulations_path::String)
        #CSV.write(simulations_path, simulations)
        print(eltype.(eachcol(simulations)))
        # open(simulations_path, "w") do f
        #     write(f, join(names(simulations), "\t") * "\n") # print header
        #     for row in 1:size(simulations)[1]
        #         line = string(simulations[row,1])
        #         for col in 2:size(simulations)[2]
        #             if typeof(simulations[row, col]) == String 
        #                 c = simulations[row, col] 
        #             else 
        #                 c = string(simulations[row, col])
        #             end
        #             line = line * "\t" * c  # merge all cells of one row
        #         end
        #         write(f, "$line\n") # print df line by line
        #     end
        # end
    end

end
