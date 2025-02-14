module Simulation_helpers

export create_simulations, create_benchmarked_simulations

using OrdinaryDiffEq
using DataFrames
using CSV
using TimerOutputs
using ModelingToolkit
using Plots
using Term

include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR, MODEL_PATH
import .Utils.Definitions: tree_definitions
import .Utils.Benchmarking: save_times_as_csv
import .Utils.Options: model_types

include("../models/julia_from_jgraph.jl")
import .Julia_from_jgraph: get_ODE_components

include(MODEL_PATH)
import .Julia_models: jf_dxdt!

# Already specified in utils.jl
trees::tree_definitions = tree_definitions()

function ODE_solver(x0, sim_options, parameters_value, sol_options)::SciMLBase.ODESolution

    prob = ODEProblem(jf_dxdt!, x0, sim_options.tspan, parameters_value)

    sol = solve(
        prob,
        sol_options.solver,
        saveat = sim_options.tpoints,
        dense = false,
        reltol = sol_options.relative_tolerance,
        abstol = sol_options.absolute_tolerance,
    )

    return sol

end

function create_simulations(g_options, sim_options, sol_options)

    graph_id::String = ""
    file_name::String = ""

    for n_node ∈ g_options.n_nodes, tree_id ∈ g_options.tree_ids
        graph_id = "$(tree_id)_$(n_node)"
        for vessel_tree ∈ trees.vascular_trees[tree_id]
            x0, p = get_ODE_components(tree_id, n_node, vessel_tree)
            simulations = ODE_solver(x0, sim_options, p, sol_options)
            if sim_options.save_simulations
                file_name = "$(graph_id)_$(vessel_tree)_$(model_type)"
                save_simulations_to_csv(
                    simulations = simulations,
                    column_names = ["$(edge[1]),$(edge[2])" for edge in p[3]],
                    simulations_path = joinpath(
                        JULIA_RESULTS_DIR,
                        tree_id,
                        graph_id,
                        "simulations",
                        "$(file_name).csv",
                    ),
                )
            end
        end
    end
end

function create_benchmarked_simulations(g_options, sim_options, bench_options, sol_options)
    to = TimerOutput()
    graph_id::String = ""

    for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
        graph_id = "$(tree_id)_$(n_node)"

        for vessel_tree ∈ trees.vascular_trees[tree_id]
            x0, p = get_ODE_components(tree_id, n_node, vessel_tree)
            printstyled("   ODEs solving started   \n"; color = 9)
            @timeit to "$(graph_id)_$(vessel_tree)" begin
                for ki in range(1, bench_options.n_iterations, step = 1)
                    @timeit to "$(graph_id)_$ki" ODE_solver(x0, sim_options, p, sol_options)
                end
            end
            (bench_options.save_running_times) && (n_species::Int64 = length(x0))
            show(to, sortby = :firstexec)
            (bench_options.save_running_times) && (save_times_as_csv(
                times = to,
                n_node = n_node,
                tree_id = "$(tree_id)_$(vessel_tree)",
                n_species = n_species,
                model_type = model_type,
                solver_name = sol_options.solver_name,
                n_term = Int32(n_species / 3),
            ))
            reset_timer!(to::TimerOutput)
        end
    end
end

function save_simulations_to_csv(;
    simulations::SciMLBase.ODESolution,
    column_names::Vector{String},
    simulations_path::String,
)
    # Convert solution to DataFrame
    df = DataFrame(simulations)

    # Write DataFrame to CSV
    header = vcat(["time"], column_names)
    CSV.write(simulations_path, df, header = header)
end
end
