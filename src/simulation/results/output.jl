"""
    simulation_path(tree_info, graph_subsystem, flow_scaling_factor, dt, format)

Path of the ARROW file with the simulation of `graph_subsystem`. The flow scaling factor
is part of the file name only if it differs from 1.0.
"""
function simulation_path(tree_info, graph_subsystem, flow_scaling_factor, dt, format::String)
    flow_part = flow_scaling_factor == 1.0 ? "" : "_$(flow_scaling_factor)Q"
    SIM_DIR = joinpath(
        tree_info.GRAPH_DIR,
        "simulations"
    )
    mkpath(SIM_DIR)
    return joinpath(
        SIM_DIR,
        "$(graph_subsystem)$(flow_part)_simulations_$(dt)_dt.$(format)",
    )
end

"""
    mmap_path(tree_info, graph_subsystem)

Path of the temporary file that holds the solution of `graph_subsystem` during the simulation
(see [`mmap_solution`](@ref)). It is in the `simulations` folder of the graph, not in `/tmp`,
which may be kept in RAM.
"""
function mmap_path(tree_info, graph_subsystem)
    return joinpath(tree_info.GRAPH_DIR, "simulations", "$(graph_subsystem)_solution.tmp")
end

"""
    check_disk_space(tree_info, n_states_total, n_time_points, sim_options)

Check the free space of the `simulations` folder of the graph against the space needed for
the solution: the temporary files (see [`mmap_solution`](@ref)) and, if
`sim_options.save_simulations`, once more the `.arrow` files (approximately, the states ids
are not counted).

- Every time point stored (`save_abstol = save_reltol = 0`): error if the `n_time_points`
  time points do not fit.
- Adaptive saving: the number of stored time points is known only after the run. Error if
  even the minimum does not fit (first and last time point, plus one every `save_max_gap`
  steps); only a warning if the worst case (all `n_time_points`) does not fit.

The temporary files are sparse, so a full disk is not detected when they are created, but only
later, when the solution is written to them: Julia then crashes with a bus error.
"""
function check_disk_space(tree_info, n_states_total::Int, n_time_points::Int, sim_options)
    simulations_dir = mkpath(joinpath(tree_info.GRAPH_DIR, "simulations"))
    n_copies = sim_options.save_simulations ? 2 : 1
    bytes(n_points) = n_copies * n_states_total * n_points * sizeof(Float64)
    gb(b) = round(b / 1e9; digits = 1)
    available = diskstat(simulations_dir).available
    worst_case = bytes(n_time_points)
    worst_case <= available && return nothing

    save_all = sim_options.save_abstol == 0 && sim_options.save_reltol == 0
    # time points that are always stored: first, last and one every save_max_gap steps
    n_min_points = save_all ? n_time_points :
        min(n_time_points, cld(n_time_points - 1, sim_options.save_max_gap) + 1)
    needed_min = bytes(n_min_points)
    if needed_min > available
        error("Not enough disk space in $(simulations_dir): " *
              "$(gb(needed_min)) GB needed at least " *
              "($(n_states_total) states, $(n_min_points) time points), " *
              "$(gb(available)) GB available")
    end
    @warn "Disk space in $(simulations_dir) may not be enough: " *
          "$(gb(worst_case)) GB needed if all $(n_time_points) time points are stored, " *
          "$(gb(available)) GB available (at most " *
          "$(fld(available, bytes(1))) time points fit; adaptive saving stores fewer)"
    return nothing
end


"""
    mmap_solution(tree_info, graph_subsystem, n_states, n_time_points) -> Matrix{Float64}

Zero matrix `(n_states, n_time_points)` backed by a temporary file in the `simulations`
folder, so that the solution of large graphs does not have to fit in RAM.
"""
function mmap_solution(tree_info, graph_subsystem, n_states, n_time_points)
    path = mmap_path(tree_info, graph_subsystem)
    M = open(io -> Mmap.mmap(io, Matrix{Float64}, (n_states, n_time_points)), path, "w+")
    return M   # file is zero-filled on creation
end


"""
    save_solutions(solutions, species_ids, times, tree_info, dt, flow_scaling_factor)

Save the solution of every subsystem as an ARROW table with one row per species: column
`species_ids`, then one column per stored time point named `t_<time>` (`times[k]` for
column `k` of `solutions[graph_subsystem]`). With adaptive saving the time points are not
equally spaced. See [`simulation_path`](@ref) for the file names.

The columns are passed to Arrow as views (`copycols = false`), without an intermediate copy.
"""
function save_solutions(solutions, species_ids, times, tree_info, dt, flow_scaling_factor)
    # rounded, so that float noise of tmin + k*dt does not end up in the column names
    column_names = ["species_ids"; ["t_$(round(t; digits = 10))" for t in times]]
    for (graph_subsystem, solution) in pairs(solutions)
        table = DataFrame(
            AbstractVector[species_ids[graph_subsystem], eachcol(solution)...],
            column_names;
            copycols = false,
        )
        save_as_arrow(
            table,
            simulation_path(tree_info, graph_subsystem, flow_scaling_factor, dt, "arrow"),
        )
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