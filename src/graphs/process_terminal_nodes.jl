module Process_Terminal_Nodes

export process_terminal_nodes

"""
Module for collecting information about terminal nodes using all individual vessel trees (arterial, portal, etc.) files.

Input: .arrow (table) for every individual vessel tree (arterial, portal, etc.).

Output: .arrow (table) with information about terminal nodes.

Idea of this module: to get, to store, and to save the information about terminal nodes that we need for correct ODEs.
    These nodes serve as a connection point between individual vessel trees.

TODO: Optimize code
"""

using ..Utils: JULIA_RESULTS_DIR
using ..Utils.Definitions: tree_definitions, ODE_groups
using ..Processing_Helpers: selection_from_df, save_as_arrow, get_extended_vector

using DataFrames, Arrow

# Already specified in utils.jl
const trees = tree_definitions()
const groups = ODE_groups()

volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000 # [cm^3] -> [l] # calculation of terminal volume



# Main function: workflow for whole tree
function process_terminal_nodes(tree_info)
    @info "Processing terminal nodes..."

    graphs = Dict{String, DataFrame}()
    collect_graphs_information!(graphs, tree_info)

    terminal_nodes_info = prepare_terminal_nodes_information(graphs, tree_info)
    save_as_arrow(terminal_nodes_info, "T", tree_info.GRAPH_DIR)
end

function collect_graphs_information!(graphs::Dict{String, DataFrame}, tree_info)
    vascular_trees = values(tree_info.vascular_trees)
    for vascular_tree ∈ Iterators.flatten(vascular_trees)
        GRAPH_PATH::String = 
            joinpath(
                tree_info.GRAPH_DIR,
                "graphs",
                "$(vascular_tree).arrow",
            )
        graphs[vascular_tree] = load_graph(GRAPH_PATH)
    end
end

function prepare_terminal_nodes_information(graphs::Dict{String, DataFrame}, tree_info)::DataFrame

    # not a good way, but for now its okay
    n_terminals::Integer = length(selection_from_df(graphs["A"], (graphs["A"].ODE_groups .== groups.terminal, :ODE_groups)))
    n_inflows::Integer = length(tree_info.vascular_trees[:inflow_trees])
    outflow_trees::Vector{String} = tree_info.vascular_trees[:outflow_trees]
    
    flow_values::Array = zeros(n_inflows+1, n_terminals)
    flow_affiliations::Array = fill([""], n_inflows+1, n_terminals)
    x_affiliations::Array = fill("", n_inflows+1, n_terminals)
    k_inflow::Integer = 2
    
    for graph in values(graphs)
        flows_single_tree = selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :flows))
        if graph[1, :is_inflow]
            flow_values[k_inflow, :] .= flows_single_tree
            flow_affiliations[k_inflow, :] .= [[graph[1, :vascular_tree_id]] for _ in 1:n_terminals]
            x_affiliations[k_inflow, :] .= [graph[1, :vascular_tree_id] for _ in 1:n_terminals]
            k_inflow += 1
        else
            flow_values[1, :] .+= flows_single_tree
        end    
        
    end

    x_affiliations[1, :] .= "T"
    flow_affiliations[1, :] .= [outflow_trees]
    volume_terminal = volume_geometry / n_terminals

    df_x_affiliations = DataFrame(x_affiliations, :auto)
    df_flow_affiliations = DataFrame(flow_affiliations, :auto)
    df_flow_values = DataFrame(flow_values, :auto)
    terminal_nodes_info = vcat(df_x_affiliations, df_flow_values, df_flow_affiliations)
    push!(terminal_nodes_info, [volume_terminal for _ in 1:n_terminals])

    return terminal_nodes_info 
end

function load_graph(GRAPH_PATH::String)::DataFrame
    # this is a duplicate of a function from julia_from_jgraph.jl
    graph = DataFrame(Arrow.Table(GRAPH_PATH))[:, [:vascular_tree_id, :is_inflow, :terminal_edges, :preterminal_edges, :flows, :ODE_groups]]
    return graph
end

end