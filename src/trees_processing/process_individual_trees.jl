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
    get_extended_vector

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
    GRAPH_PATH::String = joinpath(GRAPH_DIR, "julia/$(vascular_tree).lg")
    EDGES_PATH::String = joinpath(GRAPH_DIR, "julia/$(vascular_tree)_edges.csv")
    NODES_PATH::String = joinpath(GRAPH_DIR, "julia/$(vascular_tree)_nodes.csv")

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
        to the dataframe from which it was found out: to dataframe with edges info.
    """
    # add to the dataframe columns with species ids, flow ids, volume ids
    # we need them to be able to then save simulations with informative columns,
    # and to be able to find out for each flow and volume value to what edge they 
    # belong
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
    # add columns which will indicate whether the edge is preterminal, terminal
    # or first (start) one. We need this info to right correct equations, but
    # at this stage of graph processing, its dataframe does not have terminal edges,
    # so column :terminal is filled with "false"
    label_special_edges!(graph_structure)
    # create terminal edges and marginal edge (the one to where we add
    # intervention)
    create_special_edges!(graph_structure)
    # add column which will indicate to which ODE group (defined in utils.jl)
    # belongs each edge, storing this info in one column will fasten solving ODE
    # system
    assign_ODE_group!(graph_structure)
    # add column which will store the position (index) of each edge in the dataframe,
    # we need them to then for asch edge store indices of its preedges and postedges,
    # storing this info separately will fasten solving the ODE model
    set_index!(graph_structure)
end

function create_graph_structure(
    graph_structure::DataFrame,
    nodes_attrib::DataFrame,
    vascular_tree::String
)::DataFrame
    """
    Function that transforms columns of the dataframes with edges and nodes to separate vectors,
        prepares them to be saved as one dataframe in .arrow format, additional transforms outflow graphs.
    
    Note: .arrow format can only store tables/dataframes. To create a table/dataframe out of several vectors,
        they must have the same length. In our case these vectors have different lengths, so
        some of them must be extended with "missing" values.
    """
    # transform columns of both dataframes to separate vectors and return them in one tuple,
    # if the graph is an outflow, change source nodes to targets and vice versa, add some more additional info
    graph_info::NamedTuple =
        prepare_graph_info(graph_structure, nodes_attrib, vascular_tree)
    # number of edges graph define the number of rows in the table that will contain all info
    # about the graph that we need, so we need to store this number (see note in the describtion of this
    # function)
    df_length::Integer = length(graph_info.all_edges)
    # collect all the info in one dataframe
    # for this some vectors must be extended for them all to be the same length
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
    """Function that adds to the dataframe with edges column which indicate ODE group for each edge"""
    graph_structure.ODE_group .=
    ifelse.(
        graph_structure[:, :preterminal] .== true,
        groups.preterminal,
        ifelse.(
            graph_structure[:, :terminal] .== true,
            groups.terminal,
            ifelse.(
                graph_structure[:, :source_id] .== 0,
                groups.marginal,
                groups.other,
            ),
        ),
    )
end

function set_index!(graph_structure::DataFrame)
    """Function that adds to the dataframe with edges column which indicate position of every edge in the dataframe"""
    graph_structure.index = 1:nrow(graph_structure)
end

function prepare_graph_info(
    graph_structure::DataFrame,
    nodes_attrib::DataFrame,
    vascular_tree::String,
)::NamedTuple
    """
    Function that transforms columns of both dataframes to separate vectors, collect more info, transform outflow graph
        and returns this vectors in one tuple.
    """
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
    """
    Function that finds and returns for every edge in the graph its preedges (predecessors) and its postedges (successors).

    We need them to write ODEs correctly.
    """
    # preallocating vectors
    pre_elements = [Int64[] for _ in eachindex(edges)]
    post_elements = [Int64[] for _ in eachindex(edges)]
    # take from the dataframe with edges info only columns we need to find predecessors and successors for each edge
    df = @view graph_structure[:, [:source_id, :target_id, :index]]
    # iterate over all edges and get indices of the predecessors and successors for each edge
    @inbounds for (ke, element) in enumerate(edges)
        # unpack edge on its elements
        source_id, target_id = element
        # predecessors for the edge are edges which target ids are the same as source id of the iterated edge 
        pre_elements[ke] = selection_from_df(
            df,
            (condition.(df.target_id, source_id, ke, df.index), :index),
        )
        # successors for the edge are edges which source ids are the same as target id of the iterated edge 
        post_elements[ke] = selection_from_df(
            df,
            (condition.(df.source_id, target_id, ke, df.index), :index),
        )
    end
    return pre_elements, post_elements
end

function condition(
    column_to_look_in::I, 
    node_id_to_look::I, 
    element_index::I, 
    df_index::I
) where {I<:Integer}
    """
    Function that contains condition for the edge to be predecessor or successors of the current 
        iterated edge.
        
    :column_to_look_in - where should we search for the passed node id 
    :node_id_to_look - what node we should look for (this is a sourcs id or target id from the
        current iterated edge)
    :element_index - index (position) of the iterated edge in the vector of all edges
    :df_index - indices (positions) of the edges that we check whether they are a predecessors or 
        successors of the current iterated edge      

    Example for the predecessors:
        1. First part of the condition: check ids from the target id column in dataframe and find
            those who have the same ids as the source id of the current iterated edge
        2. Second part: becase we check whole dataframe, including terminal nodes (their source ids 
            are equal to their target ids), if we must make sure that we don't mark the current 
            iterated edge in the dataframe as predecessor or successor of itself. So,
            index of the current iterated edge must not be equal to the index of the edge we are
            checking to be predecessor or successor of it.    
    """
    return column_to_look_in == node_id_to_look && element_index != df_index
end
end
