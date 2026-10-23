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