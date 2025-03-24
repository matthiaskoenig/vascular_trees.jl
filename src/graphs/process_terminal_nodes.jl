module Process_Terminal_Nodes

export process_terminal_nodes

"""
Module for collecting information about terminal nodes using all individual vessel trees (arterial, portal, etc.) files.

Input: .arrow (table) for every individual vessel tree (arterial, portal, etc.).

Output: .arrow (table) with information about terminal nodes.

Idea of this module: to get, to store, and to save the information about terminal nodes that we need for correct ODEs.
    These nodes serve as a connection point between individual vessel trees.

TODO: Flow affiliations are wrong, DataFrame columns contain data of different types, Optimize code
"""

using ..Utils: JULIA_RESULTS_DIR
using ..Utils.Definitions: ODE_groups
using ..Processing_Helpers: selection_from_df, save_as_arrow, get_extended_vector

using DataFrames, Arrow

# Already specified in utils.jl
const groups = ODE_groups()

volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000 # [cm^3] -> [l] # calculation of terminal volume


# Main function: workflow for whole tree
function process_terminal_nodes(tree_info)
    @info "Processing terminal nodes..."

    graphs = Dict{String,DataFrame}()
    collect_graphs_information!(graphs, tree_info)

    terminal_nodes_info = prepare_terminal_nodes_information(graphs, tree_info)
    save_as_arrow(terminal_nodes_info, "T", tree_info.GRAPH_DIR)
end

function collect_graphs_information!(graphs::Dict{String,DataFrame}, tree_info)
    for vascular_tree ∈ tree_info.vascular_trees
        GRAPH_PATH::String =
            joinpath(tree_info.GRAPH_DIR, "graphs", "$(vascular_tree).arrow")
        graphs[vascular_tree] = load_graph(GRAPH_PATH)
    end
end

function prepare_terminal_nodes_information(
    graphs::Dict{String,DataFrame},
    tree_info,
)::DataFrame

    # not a good way, but for now its okay
    n_terminals::Integer = length(
        selection_from_df(
            graphs["V"],
            (graphs["V"].ODE_groups .== groups.terminal, :ODE_groups),
        ),
    )
    n_inflow::Integer = length(tree_info.tree_components[:inflow_trees])

    # preallocating vectors
    flow_values::Array = zeros(n_inflow + 1, n_terminals)
    flow_affiliations::Array{String} = fill("", n_inflow + 1, n_terminals)
    x_affiliations::Array{String} = fill("", n_inflow + 1, n_terminals)

    # where the inflow concentrations should be stored
    # row number
    k_inflow::Integer = 2

    for graph in values(graphs)
        if graph[1, :is_inflow]
            flow_values[k_inflow, :] .=
                selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :flows))
            flow_affiliations[k_inflow, :] .= selection_from_df(
                graph,
                (graph.ODE_groups .== groups.preterminal, :flow_ids),
            ) 
            x_affiliations[k_inflow, :] .= selection_from_df(
                graph,
                (graph.ODE_groups .== groups.preterminal, :species_ids),
            )
            k_inflow += 1
        else
            flow_values[1, :] .+=
                selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :flows))
        end
    end

    x_affiliations[1, :] .= ["T_$n_terminal" for n_terminal = 1:n_terminals]
    flow_affiliations[1, :] .= "Outflow"# outflow_trees
    volume_terminal = volume_geometry / n_terminals

    df_x_affiliations = DataFrame(x_affiliations, :auto)
    df_flow_affiliations = DataFrame(flow_affiliations, :auto)
    df_flow_values = DataFrame(flow_values, :auto)
    terminal_nodes_info = vcat(df_x_affiliations, df_flow_values, df_flow_affiliations)
    push!(terminal_nodes_info, [volume_terminal for _ = 1:n_terminals])

    return terminal_nodes_info
end

function load_graph(GRAPH_PATH::String)::DataFrame
    # this is a duplicate of a function from julia_from_jgraph.jl
    graph = DataFrame(Arrow.Table(GRAPH_PATH))[
        :,
        [
            :vascular_tree_id,
            :is_inflow,
            :terminal_edges,
            :preterminal_edges,
            :flows,
            :species_ids,
            :flow_ids,
            :ODE_groups,
        ],
    ]
    return graph
end

end
