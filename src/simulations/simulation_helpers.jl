module Simulation_Helpers

export create_simulations, create_benchmarked_simulations

using OrdinaryDiffEq, DataFrames, CSV, Sundials
using TimerOutputs, Plots

include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR, MODEL_PATH 
import .Utils.Definitions: tree_definitions, ODE_groups
import .Utils.Benchmarking: save_times_as_csv

# using ..Utils: JULIA_RESULTS_DIR, MODEL_PATH
# using ..Utils.Definitions: tree_definitions, ODE_groups
# using ..Utils.Benchmarking: save_times_as_csv

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
        vascular_trees = values(trees.vascular_trees[tree_id])
        inflow_trees = trees.vascular_trees[tree_id][:inflow_trees]
        outflow_trees = trees.vascular_trees[tree_id][:outflow_trees]

        x0 = Dict(vascular_tree => [0.0] for vascular_tree in Iterators.flatten(vascular_trees))
        p = Dict{String, Tuple{Bool, Vector{Float64}, Vector{Float64}, Vector{Int16}, Vector{Vector{Int64}}, Vector{Vector{Int64}}}}()
        synch_idx = Dict{String, Vector{Int64}}()

        # preparation step - get starting initial values and get parameters for all trees
        for vascular_tree ∈ Iterators.flatten(vascular_trees)
            x0[vascular_tree], p[vascular_tree] = get_ODE_components(tree_id, n_node, vascular_tree)
            if vascular_tree in trees.inflow_trees
                synch_idx[vascular_tree] = get_indices(p[vascular_tree][4], [groups.preterminal])
            else
                synch_idx[vascular_tree] = get_indices(p[vascular_tree][4], [groups.terminal])
            end
        end
        
        x0_terminal, p_terminal = get_ODE_components(tree_id, n_node, "T")
        
        t_left = tmin
        t_right = sdt
        while t_right <= tmax
            tspan .= [t_left, t_right]
            for vascular_tree ∈ Iterators.flatten(vascular_trees)
                solution = solve_ODE(x0[vascular_tree], p[vascular_tree], tspan, sol_options)
                x0[vascular_tree] = solution[end]
            end
            solution_terminal = solve_ODE(x0_terminal, p_terminal, tspan, sol_options)

            for inflow_tree in inflow_trees
                for (ki, inflow) in enumerate(p_terminal[2])
                    if inflow == inflow_tree
                        x0_terminal[:, ki] = x0[inflow_tree][synch_idx[inflow_tree]]
                    end
                end
            end
            for outflow_tree in outflow_trees
                x0[outflow_tree][synch_idx[outflow_tree]] = solution_terminal[end][:, 1]
            end

            t_left += sdt
            t_right += sdt
        end
    end
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
