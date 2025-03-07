module Processing_Helpers
"""
Helper functions used in workflow of julia graph processing.

Note on read_edges_attributes function:
Why id in this df is a target id? In Julia graph all attributes belongs to nodes.
    That is why while writing file edges_df in julia iteration goes through nnodes, but not nedges.
    Also at this step we don't care about whether graph is in inflow or an outflow,
    because the direction of the flow for every graph goes from the "highest" to "lowest" nodes
    (change of direction happens fare below this function).
    So, attributes of a node belong to the edge to this node.
    For example, 1 -> 2 -> 3. Node 1 has flow value, radius, but does not have length and pressure drop,
    because it is a start node (this is by default, check algorithm).
    Node 2 has flow value, radius, length and pressure drop ---> these values characterise edge (1, 2).
"""

using CSV, DataFrames, DataFramesMeta, InteractiveUtils, Parameters, Revise, Arrow

export read_edges,
    read_nodes_attributes,
    label_special_edges!,
    create_special_edges!,
    selection_from_df,
    create_tuples_from_dfrows,
    save_as_arrow

# calculation of terminal volume
volume_geometry = (0.100 * 0.100 * 0.10) / 1000 # [cm^3] -> [l]

function read_edges(GRAPH_PATH::String, EDGES_PATH::String)::DataFrame
    # read and prepare lg file with information about edges
    graph_structure::DataFrame = read_graph_skeleton(GRAPH_PATH)
    # read and prepare csv file with edges attributes
    edges_attrib::DataFrame = read_edges_attributes(EDGES_PATH)
    # join dfs with graph structure and edges attributes
    leftjoin!(graph_structure, edges_attrib, on = :target_id)
    disallowmissing!(graph_structure)
    
    return graph_structure
end

function read_graph_skeleton(GRAPH_PATH::String)::DataFrame
    graph_skeleton::DataFrame = CSV.read(GRAPH_PATH, DataFrame)
    select!(graph_skeleton, Not(names(graph_skeleton, Missing)))
    rename!(graph_skeleton, [:source_id, :target_id])

    return graph_skeleton
end

function read_edges_attributes(EDGES_PATH::String)::DataFrame
    edges_attrib::DataFrame = CSV.read(EDGES_PATH, DataFrame)
    @chain edges_attrib begin
        @rename! begin
            :target_id = :edge_idx
            :leaf = :leaf_count
            :radius = :radius_in_mm
            :flows = :flow_in_mm3_per_s # units are changed below
            :length = :length_in_mm
            :pressure_drop = :pressure_drop_in_kg_per_mm_s2
        end
        @transform! begin
            :flows = (:flows ./ 1000000 * 60) # Change units of the flow [mm3/s --> L/min]
            :volumes = π .* :radius .^ 2 .* :length ./ 1000000 # [mm3 --> L] 
        end
    end

    return edges_attrib
end

function read_nodes_attributes(NODES_PATH::String)::DataFrame
    nodes_attrib::DataFrame = CSV.read(NODES_PATH, DataFrame)
    @chain nodes_attrib begin
        @rename! begin
            :ids = :node_idx
            :x = :x_coord_in_mm
            :y = :y_coord_in_mm
            :z = :z_coord_in_mm
        end
    end

    return nodes_attrib
end

#=================================================================================================================================#
function label_special_edges!(graph_structure)
    label_preterminal_edges!(graph_structure)
    label_start_edges!(graph_structure)
    label_terminal_edges!(graph_structure)
end

function label_preterminal_edges!(graph_structure)
    graph_structure[!, :preterminal] =
        .!in.(graph_structure.target_id, [Set(graph_structure.source_id)])
end

function label_start_edges!(graph_structure)
    graph_structure[!, :start] =
        .!in.(graph_structure.source_id, [Set(graph_structure.target_id)])
end

function label_terminal_edges!(graph_structure)
    graph_structure[!, :terminal] = [false for _ ∈ 1:nrow(graph_structure)]
end

#=================================================================================================================================#
function create_special_edges!(graph_structure)
    create_terminal_edges!(graph_structure)
    create_marginal_edge!(graph_structure)
end

function create_terminal_edges!(graph_structure)
    # adding self edges for terminal nodes
    terminal_nodes_ids = selection_from_df(graph_structure, (graph_structure.preterminal .== true, :target_id))
    terminal_nodes_info = collect_terminal_edges_info(terminal_nodes_ids)
    append!(graph_structure, terminal_nodes_info)
end

function create_marginal_edge!(graph_structure)
    # adding marginal edge (for input)
    start_node_id = (selection_from_df(graph_structure, (graph_structure.start .== true, :source_id)))[1]
    push!(
        graph_structure,
        [
            0,
            start_node_id,
            0,
            0.0,
            0.0,
            0.0,
            0.0,
            graph_structure[graph_structure.start .== true, :volumes][1], # volume equal to the volume of start edge
            "Marginal",
            "",
            "",
            false,
            false,
            false,
        ],
    )
end

function collect_terminal_edges_info(
    terminal_node_ids,
)::DataFrame
    n_terminals = length(terminal_node_ids)
    volume_terminal = volume_geometry / n_terminals
    terminal_edges_info = DataFrame(
        :source_id => terminal_node_ids,
        :target_id => terminal_node_ids,
        :leaf .=> 0,
        ([:radius, :flows, :length, :pressure_drop] .=> 0.0)...,
        :volumes .=> volume_terminal,
        ([:element_ids, :flow_ids, :volume_ids] .=> ["CT", "QT", "VT"])...,
        ([:preterminal, :start] .=> false)...,
        :terminal .=> true,
    )
    return terminal_edges_info
end

#=================================================================================================================================#
function selection_from_df(df::AbstractDataFrame, conditions::Tuple{Union{Colon, BitVector}, Vector{Symbol}})::SubDataFrame
    return @view df[conditions...]
end

function selection_from_df(df::AbstractDataFrame, conditions::Tuple{Union{Colon, BitVector}, Symbol})::SubArray
    return @view df[conditions...]
end

function create_tuples_from_dfrows(df::AbstractDataFrame)
    return Tuple.(Tables.namedtupleiterator(df))
end

#=================================================================================================================================#
function save_as_arrow(
    graph::NamedTuple,
    tree_id::String,
    graph_id::String,
    vessel_tree::String,
    JULIA_RESULTS_DIR::String,
    folder::String
)
    ARROW_DIR::String = normpath(
        joinpath(@__FILE__, "../../..", JULIA_RESULTS_DIR, tree_id, graph_id, folder),
    )
    Arrow.write(joinpath(ARROW_DIR, "$(vessel_tree).arrow"), graph)
    
end

end
