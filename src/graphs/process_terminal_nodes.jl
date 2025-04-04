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

using DataFrames, Arrow, ProgressMeter

# Already specified in utils.jl
const groups = ODE_groups()
# volume of liver in which tree was grown
volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000 # [cm^3] -> [l] 

function process_terminal_nodes(tree_info)
    """Main function: workflow for getting and saving info about terminal nodes for the tree"""
    @info "Processing terminal nodes..."

    # initialize dictionary to store dataframes with all graphs info (container)
    graphs = Dict{String,DataFrame}()
    # collect all graphs, beloning to this tree, info 
    collect_graphs_information!(graphs, tree_info)

    terminal_nodes_info = prepare_terminal_nodes_information(graphs, tree_info)
    save_as_arrow(terminal_nodes_info, "T", tree_info.GRAPH_DIR)
end

function collect_graphs_information!(graphs::Dict{String,DataFrame}, tree_info)
    # iterate over all graphs (tree) that belongs to this tree
    for vascular_tree ∈ tree_info.vascular_trees
        # get path to the .arrow file of this currently iterated graph
        GRAPH_PATH::String =
            joinpath(tree_info.GRAPH_DIR, "graphs", "$(vascular_tree).arrow")
        # load .arrow file as dataframe and store it in the dictionary with corresponding key
        graphs[vascular_tree] = load_graph(GRAPH_PATH)
    end
end

function prepare_terminal_nodes_information(
    graphs::Dict{String,DataFrame},
    tree_info,
)::DataFrame

    @info "Preparation steps"
    # next lines are for being able to allocate vectors below
    # get the number of terminal nodes
    # not a good way, but for now its okay
    n_terminals::Integer = length(
        selection_from_df(
            graphs["V"],
            (graphs["V"].ODE_groups .== groups.terminal, :ODE_groups),
        ),
    )
    # get number of inflows for this tree
    n_inflow::Integer = length(tree_info.tree_components[:inflow_trees])
    # get number of all graphs (vascular trees)
    n_vascular_trees::Integer = length(tree_info.tree_components[:inflow_trees]) + length(tree_info.tree_components[:outflow_trees])

    # === preallocating arrays which will store values/strings that are needed to write correct ODEs for terminal part
    # or to be able to understand to which element of the model values belongs
    flow_values::Array = zeros(n_inflow + 1, n_terminals)
    flow_ids::Array{String} = fill("", n_inflow + 1, n_terminals)
    species_ids::Array{String} = fill("", n_inflow + 1, n_terminals)
    species_coordinates::Array = fill((0.0, 0.0, 0.0), n_vascular_trees, n_terminals)
    # to what tree part belong every species
    species_affiliations::Array = fill("", n_vascular_trees, n_terminals)

    # where the inflow species, concentrations should be stored (row number in an array)
    # inflows start from the row two
    k_inflow::Integer = 2

    p = Progress(length(values(graphs)); dt=0.5, color=:magenta)
    for (kg, graph) in enumerate(values(graphs))
        if graph[1, :is_inflow]
            equality_edge_to_node = 2
            flow_values[k_inflow, :] .=
                selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :flows))
            flow_ids[k_inflow, :] .= selection_from_df(
                graph,
                (graph.ODE_groups .== groups.preterminal, :flow_ids),
            ) 
            species_ids[k_inflow, :] .= selection_from_df(
                graph,
                (graph.ODE_groups .== groups.preterminal, :species_ids),
            )
            k_inflow += 1
        else
            equality_edge_to_node = 1
            flow_values[1, :] .+=
                selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :flows))
        end

        preterminal_edges = selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :all_edges))
        node_ids = selection_from_df(graph, (:, :nodes_ids))
        node_coordinates = selection_from_df(graph, (:, :nodes_coordinates))
        for (ke, edge_id) in enumerate(preterminal_edges)
            for (kn, node_id) in enumerate(skipmissing(node_ids))
                if node_id == edge_id[equality_edge_to_node]
                    species_coordinates[kg, ke] = node_coordinates[kn]
                    species_affiliations[kg, ke] = graph[1, :vascular_tree_id]
                    continue
                end
            end
        end
        next!(p)
    end

    species_ids[1, :] .= ["T_$n_terminal" for n_terminal = 1:n_terminals]
    flow_ids[1, :] .= "Outflow" # outflow_trees
    # calculation of terminal volume
    volume_terminal = volume_geometry / n_terminals

    df_species_ids = DataFrame(species_ids, :auto)
    df_species_coordinates = DataFrame(species_coordinates, :auto)
    df_species_affiliations = DataFrame(species_affiliations, :auto)
    df_flow_ids = DataFrame(flow_ids, :auto)
    df_flow_values = DataFrame(flow_values, :auto)
    terminal_nodes_info = vcat(df_species_ids, df_species_coordinates, df_species_affiliations, df_flow_ids, df_flow_values)
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
            :nodes_ids,
            :nodes_coordinates,
            :all_edges,
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
