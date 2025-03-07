module Simulation_Helpers

export create_simulations, create_benchmarked_simulations

using OrdinaryDiffEq, DataFrames, CSV, Sundials
using TimerOutputs, Plots

using ..Utils: JULIA_RESULTS_DIR, MODEL_PATH
using ..Utils.Definitions: tree_definitions, ODE_groups
using ..Utils.Benchmarking: save_times_as_csv

include("../models/julia_from_jgraph.jl")
import .Julia_from_jgraph: get_ODE_components

include("../../" * MODEL_PATH)
import .Julia_models: jf_dxdt!

# Already specified in utils.jl
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
        # save_idxs = idxs_tosave,
        #progress = true
    )

    #print(solution.destats)
    # print(solution)
    display(plot(solution))

    return solution

end

function create_simulations(g_options, sim_options, sol_options)
    tmin = sim_options.tspan[1]
    tmax = sim_options.tspan[2]
    sdt = sim_options.sdt


    for n_node ∈ g_options.n_nodes, tree_id ∈ g_options.tree_ids
        graph_id = "$(tree_id)_$(n_node)"
        vessel_trees = values(trees.vascular_trees[tree_id])
        x0 = Dict(vessel_tree => [0.0] for vessel_tree in Iterators.flatten(vessel_trees))
        p = Dict{String, Tuple{Bool, Vector{Float64}, Vector{Float64}, Vector{Int16}, Vector{Vector{Int64}}, Vector{Vector{Int64}}}}()
        synch_idx = Dict{String, Vector{Int64}}()

        # preparation step - get starting initial values and get parameters for all trees
        for vessel_tree ∈ Iterators.flatten(vessel_trees)
            x0[vessel_tree], p[vessel_tree] = get_ODE_components(tree_id, n_node, vessel_tree)
            synch_idx[vessel_tree] = get_indices(p[vessel_tree][4], [groups.terminal])
            synch_idx["T_$vessel_tree"] = get_indices(p["T"][2], ["$vessel_tree"])
        end
        x0["T"], p["T"] = get_ODE_components(tree_id, n_node, "T")
        
        t_left = tmin
        t_right = sdt
        while t_right <= tmax
            for vessel_tree ∈ Iterators.flatten(vessel_trees)
                tspan .= [t_left, t_right]
                solution = solve_single_tree(x0[vessel_tree], p[vessel_tree], tspan, sol_options)
                x0[vessel_tree] = solution[end]
            end
            t_left += sdt
            t_right += sdt
            x0["V"][synch_idx["A"]] .= x0["A"][synch_idx["A"]]
        end
    end
end

function solve_single_tree(x0, p, tspan, sol_options)

    idxs_tosave = get_indices(p[4], [groups.marginal, groups.terminal])
    simulations = solve_ODE(x0, p, tspan, sol_options)

    return simulations
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
