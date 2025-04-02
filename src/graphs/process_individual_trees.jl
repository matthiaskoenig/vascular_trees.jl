module Process_Individual_Trees

export process_julia_graph

"""
Module for processing individual vessel trees (arterial, portal, etc.).

Idea of this module: to get, to store, and to save the information about individual 
    vessel trees that we need for correct ODEs  and for creation correct terminal
    nodes file.

Input: .lg and .csv files for every individual vessel tree (arterial, portal, etc.).
        Each vessel tree (graph) has three files: 
        1 - .lg (basic structure of the graph - list of edges), 
        2 - .csv (edges and nodes info).

Output: .arrow (table) for every individual vessel tree (arterial, portal, etc.): 1 file
        for 1 vessel tree.

TODO: Ids of outflows are wrong (they are not source.id_target.id, but target.id_source.id), Optimize code
"""

using ..Utils.Definitions: flow_directions, ODE_groups
using ..Processing_Helpers:
    read_edges,
    read_nodes_attributes,
    label_special_edges!,
    create_special_edges!,
    create_tuples_from_dfrows,
    selection_from_df,
    save_as_arrow,
    get_extended_vector,
    get_path_to_file

using DataFrames, InteractiveUtils

# Already specified in utils.jl
const groups::ODE_groups = ODE_groups()
const flow_direction = flow_directions()

function process_julia_graph(tree_info)
    """Main function: workflow for whole tree"""

    @info "Processing individual trees..."

    # iterate over all parts of the tree, process them and save info
    for vascular_tree ∈ tree_info.vascular_trees
        process_individual_tree(tree_info.GRAPH_DIR, vascular_tree)
    end
end

function process_individual_tree(GRAPH_DIR::String, vascular_tree::String)
    """Workflow for individual vessel tree"""

    # initialize paths for every file that contain info about this graph
    GRAPH_PATH, EDGES_PATH, NODES_PATH = paths_initialization(GRAPH_DIR, vascular_tree)
    # load info about the graph from all three files (graph basic structure and edges info
    # are merged in one dataframe)
    graph_structure, nodes_attrib = read_graph(GRAPH_PATH, EDGES_PATH, NODES_PATH)
    # figure out additional info that we need for constructing correct ODE system and add 
    # it to the dataframe that contains edges info
    add_graph_characteristics!(graph_structure, vascular_tree)
    # collect all info from two dataframes and store it in one structure 
    graph = create_graph_structure(graph_structure, nodes_attrib, vascular_tree)
    # save structure with graph info as .arrow file
    save_as_arrow(graph, vascular_tree, GRAPH_DIR)
end

#=================================================================================================================================#
function paths_initialization(
    GRAPH_DIR::String,
    vascular_tree::String,
)::Tuple{String,String,String}
    """Function that initialises paths to graph files"""

    GRAPH_PATH::String = get_path_to_file(GRAPH_DIR, "julia", "$(vascular_tree).lg")
    EDGES_PATH::String = get_path_to_file(GRAPH_DIR, "julia", "$(vascular_tree)_edges.csv")
    NODES_PATH::String = get_path_to_file(GRAPH_DIR, "julia", "$(vascular_tree)_nodes.csv")
    return GRAPH_PATH, EDGES_PATH, NODES_PATH
end

function read_graph(
    GRAPH_PATH::String,
    EDGES_PATH::String,
    NODES_PATH::String,
)::Tuple{DataFrame,DataFrame}
    """
    Function that loads all graph files as dataframes and makes first transformations.
    """

    # read, prepare, collect first info and join files with information about edges and 
    # their attributes
    graph_structure::DataFrame = read_edges(GRAPH_PATH, EDGES_PATH)
    # read and prepare csv file with nodes attributes
    nodes_attrib::DataFrame = read_nodes_attributes(NODES_PATH)
    return graph_structure, nodes_attrib
end

function add_graph_characteristics!(graph_structure::DataFrame, vascular_tree::String)
    """
    Function that collects additional info that we need for the right ODE model and adds it
        to the dataframe from which it was found out.
    """
    transform!(
        graph_structure,
        [:source_id, :target_id] =>
            ByRow(
                (source_id, target_id) -> (
                    "$(vascular_tree)_$(source_id)_$(target_id)",
                    "Q_$(source_id)_$(target_id)",
                    "V_$(source_id)_$(target_id)",
                ),
            ) => [:species_ids, :flow_ids, :volume_ids],
    )
    label_special_edges!(graph_structure)
    create_special_edges!(graph_structure)
    assign_ODE_group!(graph_structure)
    set_index!(graph_structure)
end

function create_graph_structure(
    graph_structure::DataFrame,
    nodes_attrib::DataFrame,
    vascular_tree::String
)::DataFrame
    graph_info::NamedTuple =
        prepare_graph_info(graph_structure, nodes_attrib, vascular_tree)
    df_length::Integer = length(graph_info.all_edges)
    graph = DataFrame(
        vascular_tree_id = get_extended_vector(vascular_tree, df_length), #[vascular_tree; fill(missing, df_length-length(vascular_tree))], # vascular_tree_id
        is_inflow = get_extended_vector(graph_info.is_inflow, df_length), #[graph_info.is_inflow; fill(missing, df_length-length(graph_info.is_inflow))], # is_inflow
        nodes_ids = get_extended_vector(nodes_attrib.ids, df_length),
        nodes_coordinates = get_extended_vector(graph_info.nodes_coordinates, df_length),
        all_edges = graph_info.all_edges, # 
        terminal_edges = get_extended_vector(graph_info.terminal_edges, df_length), # terminal edges
        start_edge = get_extended_vector(graph_info.start_edge, df_length), # start edge
        preterminal_edges = get_extended_vector(graph_info.preterminal_edges, df_length), # preterminal edges
        flows = graph_structure.flows, # flows::Vector{Float64}
        volumes = graph_structure.volumes, # volumes::Vector{Float64}
        species_ids = graph_structure.species_ids,
        flow_ids = graph_structure.flow_ids,
        volume_ids = graph_structure.volume_ids,
        ODE_groups = graph_info.ODE_groups,
        pre_elements = graph_info.pre_elements,
        post_elements = graph_info.post_elements,
    )

    return graph
end

#=================================================================================================================================#
function assign_ODE_group!(graph_structure::DataFrame)
    graph_structure.ODE_group .= groups.other
    graph_structure[:, :ODE_group] .=
        ifelse.(
            graph_structure[:, :preterminal] .== true,
            groups.preterminal,
            ifelse.(
                graph_structure[:, :terminal] .== true,
                groups.terminal,
                ifelse.(
                    graph_structure[:, :source_id] .== 0,
                    groups.marginal,
                    graph_structure.ODE_group,
                ),
            ),
        )
end

function set_index!(graph_structure::DataFrame)
    graph_structure.index = 1:nrow(graph_structure)
end

function prepare_graph_info(
    graph_structure::DataFrame,
    nodes_attrib::DataFrame,
    vascular_tree::String,
)::NamedTuple
    is_inflow::Bool = in(vascular_tree, flow_direction.inflow_trees)
    # reverse edges if the tree is an outflow
    (!is_inflow) &&
        (rename!(graph_structure, [:source_id => :target_id, :target_id => :source_id]))

    edges = selection_from_df(graph_structure, (:, [:source_id, :target_id]))
    terminals = selection_from_df(
        graph_structure,
        (graph_structure.terminal .== true, [:source_id, :target_id]),
    )
    start = selection_from_df(
        graph_structure,
        (graph_structure.start .== true, [:source_id, :target_id]),
    )
    preterminals = selection_from_df(
        graph_structure,
        (graph_structure.preterminal .== true, [:source_id, :target_id]),
    )
    ODE_groups = selection_from_df(graph_structure, (:, :ODE_group))
    nodes_coord = selection_from_df(nodes_attrib, (:, [:x, :y, :z]))

    # some information must be converted to tuples
    all_edges, terminal_edges, start_edge, preterminal_edges, nodes_coordinates =
        map(create_tuples_from_dfrows, (edges, terminals, start, preterminals, nodes_coord))

    pre_elements, post_elements = get_pre_postelements(all_edges, graph_structure)

    return (
        is_inflow = is_inflow,
        nodes_coordinates = nodes_coordinates,
        all_edges = all_edges,
        terminal_edges = terminal_edges,
        start_edge = start_edge,
        preterminal_edges = preterminal_edges,
        ODE_groups = ODE_groups,
        pre_elements = pre_elements,
        post_elements = post_elements,
    )
end

function get_pre_postelements(
    edges::Vector{Tuple{T,T}},
    graph_structure::DataFrame,
) where {T<:Integer}
    pre_elements = [Int64[] for _ in eachindex(edges)]
    post_elements = [Int64[] for _ in eachindex(edges)]
    df = @view graph_structure[:, [:source_id, :target_id, :index]]
    @inbounds for (ke, element) in enumerate(edges)
        source_id, target_id = element
        pre_elements[ke] = selection_from_df(
            df,
            (condition.(df.target_id, source_id, ke, df.index), :index),
        )
        post_elements[ke] = selection_from_df(
            df,
            (condition.(df.source_id, target_id, ke, df.index), :index),
        )
    end
    return pre_elements, post_elements
end

condition(column_to_look_in, node_id_to_look, element_index, df_index) = column_to_look_in == node_id_to_look && element_index != df_index

end
