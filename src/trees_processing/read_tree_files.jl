"""
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
module ReadTreeFiles
    using CSV, DataFrames, DataFramesMeta, Parameters, Revise

    export read_edges, read_graph_skeleton, read_edges_attributes, read_nodes_attributes

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

    function read_nodes(NODES_PATH::String)::DataFrame
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

end