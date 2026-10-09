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
module ProcessingHelpers

    using CSV, DataFrames, DataFramesMeta, InteractiveUtils, Parameters, Revise, Arrow

    export label_special_edges!, create_special_edges!

    using ..TableUtils: selection_from_df

    # calculation of terminal volume
    volume_geometry = (0.100 * 0.100 * 0.10) / 1000 # [cm^3] -> [l]

    #=================================================================================================================================#
    function label_special_edges!(graph_structure::DataFrame)
        label_preterminal_edges!(graph_structure)
        label_start_edges!(graph_structure)
        label_terminal_edges!(graph_structure)
    end

    function label_preterminal_edges!(graph_structure::DataFrame)
        graph_structure[!, :preterminal] =
            .!in.(graph_structure.target_id, [Set(graph_structure.source_id)])
    end

    function label_start_edges!(graph_structure::DataFrame)
        graph_structure[!, :start] =
            .!in.(graph_structure.source_id, [Set(graph_structure.target_id)])
    end

    function label_terminal_edges!(graph_structure::DataFrame)
        graph_structure[!, :terminal] = [false for _ ∈ 1:nrow(graph_structure)]
    end

    #=================================================================================================================================#
    function create_special_edges!(graph_structure)
        create_terminal_edges!(graph_structure)
        create_marginal_edge!(graph_structure)
    end

    function create_terminal_edges!(graph_structure)
        # adding self edges for terminal nodes
        terminal_nodes_ids = selection_from_df(
            graph_structure,
            (graph_structure.preterminal .== true, :target_id),
        )
        terminal_nodes_info = collect_terminal_edges_info(terminal_nodes_ids)
        append!(graph_structure, terminal_nodes_info)
    end

    function create_marginal_edge!(graph_structure::DataFrame)
        # adding marginal edge (for input)
        start_node_id =
            (selection_from_df(graph_structure, (graph_structure.start .== true, :source_id)))[1]
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
                graph_structure[graph_structure.start.==true, :volumes][1], # volume equal to the volume of start edge
                "Marginal",
                "",
                "",
                false,
                false,
                false,
            ],
        )
    end

    function collect_terminal_edges_info(terminal_node_ids::SubArray)::DataFrame
        n_terminals = length(terminal_node_ids)
        species_ids = ["T_$n_terminal" for n_terminal = 1:n_terminals]
        flow_ids = ["QT_$n_terminal" for n_terminal = 1:n_terminals]
        volume_ids = ["VT_$n_terminal" for n_terminal = 1:n_terminals]
        volume_terminal = volume_geometry / n_terminals
        terminal_edges_info = DataFrame(
            :source_id => terminal_node_ids,
            :target_id => terminal_node_ids,
            :leaf .=> 0,
            ([:radius, :flows, :length, :pressure_drop] .=> 0.0)...,
            :volumes .=> volume_terminal,
            :species_ids => species_ids,
            :flow_ids => flow_ids,
            :volume_ids => volume_ids,
            ([:preterminal, :start] .=> false)...,
            :terminal .=> true,
        )
        return terminal_edges_info
    end

end
