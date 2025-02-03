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

    # calculation of terminal volume
    volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000  # [cm^3] -> [l]

    function read_edges(GRAPH_PATH::String, EDGES_PATH::String)::DataFrame
        # read and prepare lg file with information about edges
        graph_structure::DataFrame = read_graph_skeleton(GRAPH_PATH)
        # read and prepare csv file with edges attributes
        edges_attrib::DataFrame = read_edges_attributes(EDGES_PATH)
        # join dfs with graph structure and edges attributes
        leftjoin!(graph_structure, edges_attrib, on = :target_id)
        adding_graph_characteristics!(graph_structure)

        return graph_structure
    end

    function read_graph_skeleton(GRAPH_PATH::String)::DataFrame
        graph_skeleton::DataFrame = DataFrame(CSV.File(GRAPH_PATH))
        select!(graph_skeleton, Not(names(graph_skeleton, Missing)))
        rename!(graph_skeleton, [:source_id, :target_id])
        transform!(graph_skeleton, names(graph_skeleton, Int64) .=> ByRow(Int32), renamecols=false)

        return graph_skeleton
    end

    function read_edges_attributes(EDGES_PATH::String)::DataFrame
        edges_attrib::DataFrame = DataFrame(CSV.File(EDGES_PATH))
        @chain edges_attrib begin
            @rename! begin
                :target_id = :edge_idx
                :leaf = :leaf_count
                :radius = :radius_in_mm
                :flow = :flow_in_mm3_per_s # units are changed below [mm3/s]
                :length = :length_in_mm 
                :pressure_drop = :pressure_drop_in_kg_per_mm_s2
            end
            @transform! begin
                :target_id = Int32.(:target_id)
                :leaf = Int32.(:leaf)
                :flows = (:flow) # Change units of the flow [mm3/s --> L/min]  / 1000000 * 60
                :volumes = π .* :radius.^2 .* :length# [mm3 --> L]   / 1000000 
            end 
        end

        return edges_attrib
    end

    function adding_graph_characteristics!(graph_structure::DataFrame)
        transform!(graph_structure, [:source_id, :target_id] => ByRow((source_id, target_id) -> ["C_$(source_id)_$(target_id)", "Q_$(source_id)_$(target_id)", "V_$(source_id)_$(target_id)"]) => [:element_ids, :flow_ids, :volume_ids])

        graph_structure[!, :preterminal] = .!in.(graph_structure.target_id, [Set(graph_structure.source_id)])
        graph_structure[!, :start] = .!in.(graph_structure.source_id, [Set(graph_structure.target_id)])
        graph_structure[!, :terminal] = [false for _ ∈ 1:nrow(graph_structure)]
        @show graph_structure

        # adding self edges for terminal nodes
        preterminal_edges::SubDataFrame = subset(graph_structure, :preterminal => x -> x .== true, view=true)
        n_terminals::Int32 = nrow(preterminal_edges)
        volume_terminal::Float64 = volume_geometry / n_terminals
        terminal_node_ids = preterminal_edges[:, :target_id]
        for terminal_node_id ∈ terminal_node_ids
            push!(graph_structure, [terminal_node_id, terminal_node_id, 0, 0.0, 0.0, 0.0, 0.0, 0.0, volume_terminal,"CT_$(terminal_node_id)", "QT_$(terminal_node_id)", "VT_$(terminal_node_id)", false, false, true])
        end

        # adding marginal edge (for input)
        start_edges::SubDataFrame = subset(graph_structure, :start => x -> x .== true, view=true)
        start_nodes_ids = start_edges[:, :source_id]
        for start_node_id ∈ start_nodes_ids
            push!(graph_structure, [0, start_node_id, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,"Marginal", "", "", false, false, false])
        end

    end

    function read_nodes_attributes(NODES_PATH::String)::DataFrame
        nodes_attrib::DataFrame = DataFrame(CSV.File(NODES_PATH))
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