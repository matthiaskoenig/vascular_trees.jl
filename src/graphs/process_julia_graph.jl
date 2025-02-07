module Process_julia_graph
include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Definitions: tree_definitions, ODE_groups
import .Utils.Options: graph_options
import .Utils.Nice_printing: print_graph_info

include("./processing_helpers.jl")
import .Processing_helpers:
    read_edges,
    read_nodes_attributes,
    label_special_edges!,
    create_special_edges!,
    create_tuples_from_dfrows,
    selection_from_df

using DataFrames, InteractiveUtils, Arrow

"""
Module for processing julia graph files from SyntheticVascularTrees.jl.
    Only works with "Rectangle_quad" and "Rectangle_trio". Otherwise "tree_definitions" must be modified in Utils.jl.

Input: .lg and .csv files for every individual vessel tree (arterial, portal, etc.)

Output: .arrow in form of the table for every individual vessel tree (arterial, portal, etc.)

Idea of this module: to get, to store, and to save all the information that we need for correct ODEs from the graph. 

TODO: Get rid of zeros in DataFrame
"""

# ============ Specify options
g_options::graph_options = graph_options(
    n_nodes = [10, 30, 50, 100, 250, 500, 750, 1000, 1250, 1500, 1750, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 200000, 300000, 400000, 500000, 1000000],  #750, 1000, 1250, 1500
    tree_ids = [
        "Rectangle_quad",
        # "Rectangle_trio",
    ],
)

# Already specified in utils.jl
trees::tree_definitions = tree_definitions()
groups::ODE_groups = ODE_groups()

function __init__()
    for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
        process_julia_graph(tree_id, n_node)
    end
end

# Main function: workflow for whole tree
function process_julia_graph(tree_id::String, n_node::Int32)
    # get graph id for correct path definition
    graph_id::String = "$(tree_id)_$(n_node)"
    # print tree info in terminal
    print_graph_info(graph_id, tree_id)
    # get directory with graph's files
    GRAPH_DIR::String = normpath(
        joinpath(@__FILE__, "../../..", JULIA_RESULTS_DIR, tree_id, graph_id, "julia"),
    )
    for vessel_tree ∈ trees.vascular_trees[tree_id]
        process_individual_tree(GRAPH_DIR, vessel_tree, tree_id, graph_id)
    end
end

# Workflow for individual vessel trees
function process_individual_tree(
    GRAPH_DIR::String,
    vessel_tree::String,
    tree_id::String,
    graph_id::String,
)
    GRAPH_PATH, EDGES_PATH, NODES_PATH = paths_initialization(GRAPH_DIR, vessel_tree)
    graph_structure, nodes_attrib = read_graph(GRAPH_PATH, EDGES_PATH, NODES_PATH)
    add_graph_characteristics!(graph_structure)
    graph = create_graph_structure(graph_structure, nodes_attrib, vessel_tree)
    save_graph(graph, tree_id, graph_id, vessel_tree)
end

#=================================================================================================================================#
function paths_initialization(
    GRAPH_DIR::String,
    vessel_tree::String,
)::Tuple{String,String,String}
    GRAPH_PATH::String = joinpath(GRAPH_DIR, "$(vessel_tree).lg")
    EDGES_PATH::String = joinpath(GRAPH_DIR, "$(vessel_tree)_edges.csv")
    NODES_PATH::String = joinpath(GRAPH_DIR, "$(vessel_tree)_nodes.csv")
    return GRAPH_PATH, EDGES_PATH, NODES_PATH
end

function read_graph(
    GRAPH_PATH::String,
    EDGES_PATH::String,
    NODES_PATH::String,
)::Tuple{DataFrame,DataFrame}
    # read, prepare and join files with information about edges and their attributes
    graph_structure::DataFrame = read_edges(GRAPH_PATH, EDGES_PATH)
    # read and prepare csv file with nodes attributes
    nodes_attrib::DataFrame = read_nodes_attributes(NODES_PATH)

    return graph_structure, nodes_attrib
end

function add_graph_characteristics!(graph_structure::DataFrame)
    transform!(
        graph_structure,
        [:source_id, :target_id] =>
            ByRow(
                (source_id, target_id) -> (
                    "C_$(source_id)_$(target_id)",
                    "Q_$(source_id)_$(target_id)",
                    "V_$(source_id)_$(target_id)",
                ),
            ) => [:element_ids, :flow_ids, :volume_ids],
    )
    label_special_edges!(graph_structure)
    create_special_edges!(graph_structure)
    assign_ODE_group!(graph_structure)
    set_index!(graph_structure)
end

function create_graph_structure(
    graph_structure::DataFrame,
    nodes_attrib::DataFrame,
    vessel_tree::String,
)::DataFrame
    graph_info::NamedTuple = prepare_graph_info(graph_structure, nodes_attrib, vessel_tree)
    df_length = maximum(length, [graph_info.all_edges, nodes_attrib.ids])
    graph = DataFrame(
        vascular_tree_id = [vessel_tree; fill(missing, df_length-length(vessel_tree))], # vascular_tree_id
        is_inflow = [graph_info.is_inflow; fill(missing, df_length-length(graph_info.is_inflow))], # is_inflow
        nodes_ids = [nodes_attrib.ids; fill(missing, df_length-length(nodes_attrib.ids))], # nodes::Vector{Int32}
        nodes_coordinates = [graph_info.nodes_coordinates; fill(missing, df_length-length(graph_info.nodes_coordinates))],
        all_edges = [graph_info.all_edges; fill(missing, df_length-length(graph_info.all_edges))], # 
        terminal_edges = [graph_info.terminal_edges; fill(missing, df_length-length(graph_info.terminal_edges))], # terminal edges
        start_edge = [graph_info.start_edge; fill(missing, df_length-length(graph_info.start_edge))], # start edge
        preterminal_edges = [graph_info.preterminal_edges; fill(missing, df_length-length(graph_info.preterminal_edges))], # preterminal edges
        flows = [graph_structure.flows; fill(missing, df_length-length(graph_structure.flows))], # flows::Vector{Float64}
        volumes = [graph_structure.volumes; fill(missing, df_length-length(graph_structure.volumes))], # volumes::Vector{Float64}
        element_ids = [graph_structure.element_ids; fill(missing, df_length-length(graph_structure.element_ids))],
        flow_ids = [graph_structure.flow_ids; fill(missing, df_length-length(graph_structure.flow_ids))],
        volume_ids = [graph_structure.volume_ids; fill(missing, df_length-length(graph_structure.volume_ids))],
        ODE_groups =[graph_info.ODE_groups; fill(missing, df_length-length(graph_info.ODE_groups))],
        pre_elements = [graph_info.pre_elements; fill(missing, df_length-length(graph_info.pre_elements))],
        post_elements = [graph_info.post_elements; fill(missing, df_length-length(graph_info.post_elements))],
    )
    return graph
end

function save_graph(
    graph::DataFrame,
    tree_id::String,
    graph_id::String,
    vessel_tree::String,
)
    ARROW_DIR::String = normpath(
        joinpath(@__FILE__, "../../..", JULIA_RESULTS_DIR, tree_id, graph_id, "graphs"),
    )
    Arrow.write(joinpath(ARROW_DIR, "$(vessel_tree).arrow"), graph)
    
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
    graph_structure.index = Int32.(1:nrow(graph_structure))
end

function prepare_graph_info(graph_structure::DataFrame, nodes_attrib::DataFrame, vessel_tree::String)::NamedTuple
    is_inflow::Bool = in(vessel_tree, trees.inflow_trees)
    # reverse edges if the tree is an outflow
    (!is_inflow) &&
        (rename!(graph_structure, [:source_id => :target_id, :target_id => :source_id]))

    edges = selection_from_df(graph_structure, (:, [:source_id, :target_id]))
    terminals = selection_from_df(graph_structure, (graph_structure.terminal .== true, [:source_id, :target_id]))
    start = selection_from_df(graph_structure, (graph_structure.start .== true, [:source_id, :target_id]))
    preterminals = selection_from_df(graph_structure, (graph_structure.preterminal .== true, [:source_id, :target_id]))
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
    edges::Vector{Tuple},
    graph_structure::DataFrame,
)::NTuple{2,Vector{Vector{Int32}}}
    pre_elements = [[] for _ in eachindex(edges)]
    post_elements = [[] for _ in eachindex(edges)]
    @inbounds for (ke, element) in enumerate(edges)
        source_id, target_id = element
        pre_elements[ke] = graph_structure[
            graph_structure.target_id.==source_id .&& ke.!=graph_structure.index,
            :index,
        ]
        post_elements[ke] = graph_structure[
            graph_structure.source_id.==target_id .&& ke.!=graph_structure.index,
            :index,
        ]
    end
    return pre_elements, post_elements
end

end
