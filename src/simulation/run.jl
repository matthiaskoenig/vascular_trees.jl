"""
    run_simulations(tree_info, sim_options, sol_options, additional_sol_options,
                    flow_scaling_factor, bench_options)

Create and run the coupled tree simulation for one graph and one flow scaling factor.

If `sim_options.benchmark` is `true`, the simulation is run `bench_options.n_iterations` times
and timed (the results are printed and, if `bench_options.save_running_times`, saved to a CSV).
Otherwise it is run once and, if `sim_options.save_simulations`, the solutions are saved to ARROW.

# Arguments
- `tree_info`: `TreeStructure` of the graph (ids of the trees, graph directory, ...).
- `sim_options`: simulation options (`tspan`, `steps`, `dt`, `save_simulations`, `benchmark`).
- `sol_options`: ODE solver options (solver, tolerances, `compute_jacobian`).
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
    trees = tree_info.vascular_trees
    n_time_points::Int = sim_options.steps + 1

    p          = Dictionary{String, VascularTreeParameters{Float64, Int}}()
    u0         = Dictionary{String, Vector{Float64}}()
    synch_idxs = Dictionary{String, Vector{Int}}()
    solutions  = Dictionary{String, Matrix{Float64}}()

    # terminal nodes are stored separately from trees, because they are stored differently
    p_terminal = get_ODE_parameters(tree_info, "T", flow_scaling_factor)
    # the state is the (n_rows, n_terminals) block as a vector (column-major order, see
    # jf_terminal!): the analytic-Jacobian callback of Sundials.jl only supports vectors
    u0_terminal = vec(get_initial_values(p_terminal.flow_values))

    for (tree_row, tree) in enumerate(trees)
        # 1. parameters of the tree
        p_tree = get_ODE_parameters(tree_info, tree, flow_scaling_factor)
        insert!(p, tree, p_tree)
        # 2. initial values and indices of the species that connect trees and terminal nodes;
        #    row tree_row of synch_indices belongs to this tree (order of tree_info.vascular_trees)
        insert!(u0, tree, get_initial_values(p_tree.ODE_groups))
        insert!(synch_idxs, tree,
            get_synchronization_indices(tree, p_tree.ODE_groups, view(p_terminal.synch_indices, tree_row, :)))
    end

    # 3. preallocated solutions, one column per (possibly) stored time point; the disk space
    #    is checked first, because a full disk is noticed only while writing (crash)
    # own keys: map(length, u0) would share the keys of u0, and the insert! below would add
    # "T" (without a value) to u0 as well
    n_states = dictionary(id => length(u0[id]) for id in keys(u0))
    insert!(n_states, "T", length(u0_terminal))
    if sim_options.benchmark
        # nothing is stored when benchmarking (see integrate!): empty matrices, no temporary files
        for (graph_subsystem, n) in pairs(n_states)
            insert!(solutions, graph_subsystem, Matrix{Float64}(undef, n, 0))
        end
    else
        stage("preallocating solutions ($(sum(n_states)) states * $n_time_points time points)") do
            check_disk_space(tree_info, sum(n_states), n_time_points, sim_options.save_simulations)
            for (graph_subsystem, n) in pairs(n_states)
                insert!(solutions, graph_subsystem, mmap_solution(tree_info, graph_subsystem, n, n_time_points))
            end
        end
    end

    try
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
                species_counts = count_species(p_terminal, p)
                # one column per field, see save_times_as_csv
                run_info = (
                    n_states           = species_counts.n_states,
                    n_distinct_species = species_counts.n_distinct,
                    solver             = sol_options.solver_name,
                    n_node             = tree_info.n_node,
                    tree_id            = tree_info.graph_id,
                    n_terminals        = size(p_terminal.flow_values, 2),
                    tolerance          = sol_options.relative_tolerance,
                    abs_tolerance      = sol_options.absolute_tolerance,
                    jacobian           = sol_options.compute_jacobian,
                    dt                 = sim_options.dt,
                    save_reltol        = sim_options.save_reltol,
                    save_abstol        = sim_options.save_abstol,
                    version            = bench_options.version,
                )
                save_times_as_csv(to, run_info)
            end
            reset_timer!(to)
        else
            # stored time points (adaptive saving, see integrate!)
            times = solve_tree!(
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
            @info "Stored $(length(times)) of $(n_time_points) time points"
            if sim_options.save_simulations
                # 4. species ids — same keys as solutions; terminal ids in linear (column-major)
                #    order, the order of the terminal state vector
                species_ids = similar(solutions, Vector{String})
                for vascular_tree ∈ tree_info.vascular_trees
                    species_ids[vascular_tree] = p[vascular_tree].species_ids
                end
                species_ids["T"] = vec(p_terminal.species_ids)
                # only the stored columns, the rest of the preallocated matrix is unused
                stored_solutions = map(solution -> view(solution, :, 1:length(times)), solutions)
                stage("saving results") do
                    save_solutions(stored_solutions, species_ids, times, tree_info, sim_options.dt, flow_scaling_factor)
                end
            end
        end
    finally
        # temporary storage of the solutions (see mmap_solution), the results are in the .arrow files
        for graph_subsystem in keys(solutions)
            rm(mmap_path(tree_info, graph_subsystem); force = true)
        end
    end
end

"""
    count_species(p_terminal, p) -> (n_states, n_distinct)

Number of species of the coupled system, counted in two ways.
- `n_states`: states integrated (and stored) by the solvers: the whole terminal block plus
  all states of the trees, i.e. the sum of the lengths of all `u0` vectors. It contains
  copies: rows 2..n of the terminal block are copies of the inflow preterminal species, and
  every tree has one terminal self-edge per terminal, a copy of the terminal compartment
  (both set by synchronization).
- `n_distinct`: physically distinct species: only the terminal compartments (row 1 of the
  terminal block) plus the tree states without terminal edges.
"""
function count_species(p_terminal, p)
    n_tree_states  = sum(p_tree -> length(p_tree.ODE_groups), p)
    n_tree_species = sum(p_tree -> count(!=(groups.terminal), p_tree.ODE_groups), p)
    n_rows, n_terminals = size(p_terminal.flow_values)
    return (n_states = n_rows * n_terminals + n_tree_states,
            n_distinct = n_terminals + n_tree_species)
end
