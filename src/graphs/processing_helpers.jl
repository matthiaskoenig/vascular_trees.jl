module Processing_helpers
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

    using CSV, DataFrames, DataFramesMeta, InteractiveUtils

    export read_edges, read_nodes_attributes

    function read_edges(GRAPH_DIR::String)::DataFrame
        # read and prepare lg file with information about edges
        graph_structure::DataFrame = read_graph_skeleton(GRAPH_DIR)
        # read and prepare csv file with edges attributes
        edges_attrib::DataFrame = read_edges_attributes(GRAPH_DIR)
        # join dfs with graph structure and edges attributes
        leftjoin!(graph_structure, edges_attrib, on = :target_id)

        return graph_structure
    end

    function read_graph_skeleton(GRAPH_DIR::String)::DataFrame
        graph_skeleton::DataFrame = DataFrame(CSV.File(joinpath(GRAPH_DIR, "A.lg")))
        select!(graph_skeleton, Not(names(graph_skeleton, Missing)))
        rename!(graph_skeleton, [:source_id, :target_id])
        transform!(graph_skeleton, names(graph_skeleton, Int64) .=> ByRow(Int32), renamecols=false)
        graph_skeleton[!, :terminal] = .!in.(graph_skeleton.target_id, [Set(graph_skeleton.source_id)])
        graph_skeleton[!, :start] = .!in.(graph_skeleton.source_id, [Set(graph_skeleton.target_id)])

        terminal_edges::SubDataFrame = subset(graph_skeleton, :terminal => x -> x .== true, view=true)
        graph_skeleton[!, :preterminal] = in.(graph_skeleton.target_id, [Set(terminal_edges.source_id)])

        return graph_skeleton
    end

    function read_edges_attributes(GRAPH_DIR::String)::DataFrame
        edges_attrib::DataFrame = DataFrame(CSV.File(joinpath(GRAPH_DIR, "A_edges.csv")))
        @chain edges_attrib begin
            @rename! begin
                :target_id = :edge_idx
                :leaf = :leaf_count
                :radius = :radius_in_mm
                :flow = :flow_in_mm3_per_s # units are changed below
                :length = :length_in_mm 
                :pressure_drop = :pressure_drop_in_kg_per_mm_s2
            end
            @transform! begin
                :target_id = Int32.(:target_id)
                :leaf = Int32.(:leaf)
                :flows = (:flow ./ 1e6 .* 60) # Change units of the flow [mm3/s --> L/min]
                :volumes = π .* :radius.^2 .* :length ./ 1e6 # [mm3 --> L]
            end 
        end

        return edges_attrib
    end

    function create_terminal_edges(graph_structure::DataFrame)
        terminal_edges::SubDataFrame = subset(graph_skeleton, :terminal => x -> x .== true, view=true)
    end

    function read_nodes_attributes(GRAPH_DIR::String)::DataFrame
        nodes_attrib::DataFrame = DataFrame(CSV.File(joinpath(GRAPH_DIR, "A_nodes.csv")))
        @chain nodes_attrib begin
            @rename! begin
                :ids = :node_idx
                :x = :x_coord_in_mm
                :y = :y_coord_in_mm
                :z = :z_coord_in_mm
            end
            @transform! :ids = Int32.(:ids)
        end

        return nodes_attrib
    end
end