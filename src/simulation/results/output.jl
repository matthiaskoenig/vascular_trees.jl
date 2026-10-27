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
    mmap_path(tree_info, graph_subsystem)

Path of the temporary file that holds the solution of `graph_subsystem` during the simulation
(see [`mmap_solution`](@ref)). It is in the `simulations` folder of the graph, not in `/tmp`,
which may be kept in RAM.
"""
function mmap_path(tree_info, graph_subsystem)
    return joinpath(tree_info.GRAPH_DIR, "simulations", "$(graph_subsystem)_solution.tmp")
end

"""
    check_disk_space(tree_info, n_species_total, n_time_points, save_results)

Throw an error if the `simulations` folder of the graph has less free space than the worst
case: every time point stored in the temporary files (see [`mmap_solution`](@ref)) and, if
`save_results`, once more in the `.arrow` files (approximately, the species ids are not counted).

The temporary files are sparse, so a full disk is not detected when they are created, but only
later, when the solution is written to them: Julia then crashes with a bus error.
"""
function check_disk_space(tree_info, n_species_total::Int, n_time_points::Int, save_results::Bool)
    simulations_dir = mkpath(joinpath(tree_info.GRAPH_DIR, "simulations"))
    n_copies = save_results ? 2 : 1
    needed = n_copies * n_species_total * n_time_points * sizeof(Float64)
    available = diskstat(simulations_dir).available
    if needed > available
        error("Not enough disk space in $(simulations_dir): " *
              "$(round(needed / 1e9; digits = 1)) GB needed in the worst case " *
              "($(n_species_total) species, $(n_time_points) time points), " *
              "$(round(available / 1e9; digits = 1)) GB available")
    end
end


"""
    mmap_solution(tree_info, graph_subsystem, n_species, n_time_points) -> Matrix{Float64}

Zero matrix `(n_species, n_time_points)` backed by a temporary file in the `simulations`
folder, so that the solution of large graphs does not have to fit in RAM.
"""
function mmap_solution(tree_info, graph_subsystem, n_species, n_time_points)
    path = mmap_path(tree_info, graph_subsystem)
    M = open(io -> Mmap.mmap(io, Matrix{Float64}, (n_species, n_time_points)), path, "w+")
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