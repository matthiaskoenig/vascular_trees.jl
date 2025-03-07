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
using ..Processing_Helpers: selection_from_df, save_as_arrow

using DataFrames, Arrow

# Already specified in utils.jl
const trees::tree_definitions = tree_definitions()
const groups::ODE_groups = ODE_groups()


# Main function: workflow for whole tree
function process_terminal_nodes(tree_id::String, n_node::Integer)
    @info "Processing terminal nodes..."
    # get graph id for correct path definition
    graph_id::String = "$(tree_id)_$(n_node)"
    vessel_trees = values(trees.vascular_trees[tree_id])

    graphs = Dict{String, DataFrame}()
    collect_graphs_information!(graphs, vessel_trees, graph_id, tree_id)
    terminal_nodes_info = prepare_terminal_nodes_information(graphs)
    save_as_arrow(terminal_nodes_info, tree_id, graph_id, "T", JULIA_RESULTS_DIR, "graphs")
end

function collect_graphs_information!(graphs::Dict{String, DataFrame}, vessel_trees::Base.ValueIterator, graph_id::String, tree_id::String)
    for vessel_tree ∈ Iterators.flatten(vessel_trees)
        GRAPH_PATH::String = normpath(
            joinpath(
                @__FILE__,
                "../../..",
                JULIA_RESULTS_DIR,
                tree_id,
                graph_id,
                "graphs/$(vessel_tree).arrow",
            ),
        )
        graphs[vessel_tree] = load_graph(GRAPH_PATH)
    end
end

function prepare_terminal_nodes_information(graphs::Dict{String, DataFrame})::NamedTuple
    flow_values::Vector{AbstractFloat} = []
    flow_affiliations::Vector{String} = []
    x_affiliations::Vector{String} = []
    for graph in values(graphs)
        flows_single_tree = selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :flows))
        push!(flow_values, flows_single_tree...)
        push!(flow_affiliations, [graph[1, :vascular_tree_id] for _ in eachindex(flows_single_tree)]...)
        if graph[1, :is_inflow] 
            x_affiliation = [graph[1, :vascular_tree_id] for _ in eachindex(flows_single_tree)]
        else
            x_affiliation = ["T" for _ in eachindex(skipmissing(graph[:, :terminal_edges]))]
        end
        push!(x_affiliations, x_affiliation...)
    end
    terminal_nodes_info = (
        x_affiliations=x_affiliations,
        flow_affiliations=flow_affiliations,
        flow_values=flow_values
    )
    return terminal_nodes_info 
end


function load_graph(GRAPH_PATH::String)::DataFrame
    graph = DataFrame(Arrow.Table(GRAPH_PATH))[:, [:vascular_tree_id, :is_inflow, :terminal_edges, :preterminal_edges, :flows, :ODE_groups]]
    return graph
end

end