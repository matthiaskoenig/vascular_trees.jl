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

# calculation of terminal volume
volume_geometry = (0.100 * 0.100 * 0.10) / 1000 # [cm^3] -> [l]


# Main function: workflow for whole tree
function process_terminal_nodes(tree_id::String, n_node::Integer)
    @info "Processing terminal nodes..."
    # get graph id for correct path definition
    graph_id::String = "$(tree_id)_$(n_node)"
    vascular_trees = values(trees.vascular_trees[tree_id])
    # get number of inflow trees
    n_inflows = length(trees.vascular_trees[tree_id][:inflow_trees])
    #get name of outflows
    outflow_trees = trees.vascular_trees[tree_id][:outflow_trees]

    graphs = Dict{String, DataFrame}()
    collect_graphs_information!(graphs, vascular_trees, graph_id, tree_id)

    terminal_nodes_info = prepare_terminal_nodes_information(graphs, n_inflows, outflow_trees)
    save_as_arrow(terminal_nodes_info, tree_id, graph_id, "T", JULIA_RESULTS_DIR, "graphs")
end

function collect_graphs_information!(graphs::Dict{String, DataFrame}, vascular_trees::Base.ValueIterator, graph_id::String, tree_id::String)
    for vascular_tree ∈ Iterators.flatten(vascular_trees)
        GRAPH_PATH::String = normpath(
            joinpath(
                @__FILE__,
                "../../..",
                JULIA_RESULTS_DIR,
                tree_id,
                graph_id,
                "graphs/$(vascular_tree).arrow",
            ),
        )
        graphs[vascular_tree] = load_graph(GRAPH_PATH)
    end
end

function prepare_terminal_nodes_information(graphs::Dict{String, DataFrame}, n_inflows::Integer, outflow_trees::Vector{String})::DataFrame

    # not a good way, but for now its okay
    n_terminals::Integer = length(selection_from_df(graphs["A"], (graphs["A"].ODE_groups .== groups.terminal, :ODE_groups)))
    
    flow_values::Array = zeros(n_inflows+1, n_terminals)
    flow_affiliations::Array = fill([""], n_inflows+1, n_terminals)
    x_affiliations::Array = fill("", n_inflows+1, n_terminals)
    k_inflow::Integer = 2
    
    for graph in values(graphs)
        flows_single_tree = selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :flows))
        if graph[1, :is_inflow]
            flow_values[k_inflow, :] = flows_single_tree
            flow_affiliations[k_inflow, :] .= [[graph[1, :vascular_tree_id]] for _ in 1:n_terminals]
            x_affiliations[k_inflow, :] = [graph[1, :vascular_tree_id] for _ in 1:n_terminals]
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
    display(terminal_nodes_info)

    return terminal_nodes_info 
end

function load_graph(GRAPH_PATH::String)::DataFrame
    # this is a duplicate of a function from julia_from_jgraph.jl
    graph = DataFrame(Arrow.Table(GRAPH_PATH))[:, [:vascular_tree_id, :is_inflow, :terminal_edges, :preterminal_edges, :flows, :ODE_groups]]
    return graph
end

end