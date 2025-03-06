module Simulation_helpers

export create_simulations, create_benchmarked_simulations

using OrdinaryDiffEq
using DataFrames
using CSV
using TimerOutputs
using ModelingToolkit
using Plots
using Term
using Sundials
# using Distributed

include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR, MODEL_PATH
import .Utils.Definitions: tree_definitions, ODE_groups
import .Utils.Benchmarking: save_times_as_csv

include("../models/julia_from_jgraph.jl")
import .Julia_from_jgraph: get_ODE_components

include("../../" * MODEL_PATH)
import .Julia_models: jf_dxdt!

# Already specified in utils.jl
const trees = tree_definitions()
const groups = ODE_groups()

function solve_ODE(x0, p, sim_options, sol_options, idxs_tosave)

    problem = ODEProblem(jf_dxdt!, x0, sim_options.tspan, p)

    solution = solve(
        problem,
        sol_options.solver,
        # saveat = sim_options.tpoints,
        dense = false,
        reltol = sol_options.relative_tolerance,
        abstol = sol_options.absolute_tolerance,
        save_idxs = idxs_tosave,
        #progress = true
    )

    #print(solution.destats)
    #display(plot(solution))

    return solution

end

function create_simulations(g_options, sim_options, sol_options)

    for n_node ∈ g_options.n_nodes, tree_id ∈ g_options.tree_ids
        graph_id = "$(tree_id)_$(n_node)"
        vessel_trees = values(trees.vascular_trees[tree_id])

        for vessel_tree ∈ Iterators.flatten(vessel_trees)
            x0, graph_p = get_ODE_components(tree_id, n_node, vessel_tree)
            p = (
                graph_p.is_inflow,
                graph_p.flows,
                graph_p.volumes,
                graph_p.ODE_groups,
                graph_p.pre_elements,
                graph_p.post_elements,
            )
            idxs_tosave = get_indices(p[4], [groups.marginal, groups.terminal])
            simulations = solve_ODE(x0, p, sim_options, sol_options, idxs_tosave)
            display(plot(simulations))

            if sim_options.save_simulations
                file_name = "$(graph_id)_$(vessel_tree)"
                save_simulations_to_csv(
                    simulations,
                    ["$(edge[1]),$(edge[2])" for edge in graph_p.all_edges], # column names
                    joinpath(
                        JULIA_RESULTS_DIR,
                        tree_id,
                        graph_id,
                        "simulations",
                        "$(file_name).csv",
                    ), # simulations_path
                )
            end
        end
    end
end

function create_simulations(g_options, sim_options, bench_options, sol_options)
    to = TimerOutput()
    graph_id::String = ""

    for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
        graph_id = "$(tree_id)_$(n_node)"
        vessel_trees = values(trees.vascular_trees[tree_id])

        for vessel_tree ∈ Iterators.flatten(vessel_trees)
            x0, graph_p = get_ODE_components(tree_id, n_node, vessel_tree)
            p = (
                graph_p.is_inflow,
                graph_p.flows,
                graph_p.volumes,
                graph_p.ODE_groups,
                graph_p.pre_elements,
                graph_p.post_elements,
            )
            idxs_tosave = get_indices(p[4], [groups.marginal, groups.terminal])
            @timeit to "$(graph_id)_$(vessel_tree)" begin
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
                    "$(tree_id)_$(vessel_tree)",
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
