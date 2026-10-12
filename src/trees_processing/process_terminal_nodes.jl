"""
Module for collecting information about terminal nodes using all individual vessel trees (arterial, portal, etc.) files.

Input: .arrow (table) for every individual vessel tree (arterial, portal, etc.).

Output: .arrow (table) with information about terminal nodes.

Idea of this module: to get, to store, and to save the information about terminal nodes that we need for correct ODEs.
    These nodes serve as a connection point between individual vessel trees.

TODO: Flow affiliations are wrong, DataFrame columns contain data of different types, Optimize code

TODO: prepare_tnodes_info assumes that every tree lists its preterminal edges in the same order,
    so that the k-th preterminal edge of each tree belongs to terminal k (column k of the output).
    This is not checked. If it does not hold, values of different terminals get mixed up
    (possible cause of the wrong flow affiliations).
"""
module ProcessTerminalNodes

    export process_terminal_nodes

    using ..Paths: JULIA_RESULTS_DIR
    using ..Definitions: ODE_groups
    using ..TableUtils: 
        selection_from_df, 
        save_as_arrow, 
        get_extended_vector,
        load_arrow

    using DataFrames,
        Arrow,
        ProgressMeter,
        OrderedCollections

    const groups = ODE_groups()
    # Volume of the liver cube in which the tree was grown (0.1 cm sides): [cm^3] -> [L].
    # It is split equally between all terminal nodes.
    volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000

    """
        process_terminal_nodes(tree_info)

    Load every individual vessel tree listed in `tree_info.vascular_trees`, collect the
    terminal node information with [`prepare_tnodes_info`](@ref), and save it to
    `<tree_info.GRAPH_DIR>/graphs/T.arrow`.

    `tree_info` must provide `GRAPH_DIR`, `vascular_trees` and `tree_components`
    (with the keys `:inflow_trees` and `:outflow_trees`).
    """
    function process_terminal_nodes(tree_info)
        @info "Processing terminal nodes..."

        # ordered, so that rows of the output follow the order of tree_info.vascular_trees
        graphs = OrderedDict{String,DataFrame}()
        # collect info from all graphs, beloning to this tree 
        collect_graphs_information!(graphs, tree_info)

        tnodes_info_df = prepare_tnodes_info(graphs, tree_info)
        save_as_arrow(tnodes_info_df, joinpath(tree_info.GRAPH_DIR, "graphs/T.arrow"))
    end

    """
        collect_graphs_information!(graphs::OrderedDict{String,DataFrame}, tree_info)

    Load the `.arrow` table of every tree in `tree_info.vascular_trees` into `graphs`,
    keyed by the tree name (e.g. `"A"`, `"V"`). Only the columns needed for the terminal
    nodes are loaded.
    """
    function collect_graphs_information!(graphs::OrderedDict{String,DataFrame}, tree_info)
        # iterate over all graphs (tree) that belongs to this tree
        for vascular_tree ∈ tree_info.vascular_trees
            # get path to the .arrow file of this currently iterated graph
            GRAPH_PATH::String =
                joinpath(tree_info.GRAPH_DIR, "graphs", "$(vascular_tree).arrow")
            # load .arrow file as dataframe and store it in the dictionary with corresponding key
            graphs[vascular_tree] = load_arrow(
                GRAPH_PATH, 
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
                ]
            )
        end
    end

    """
        prepare_tnodes_info(graphs::OrderedDict{String,DataFrame}, tree_info) -> DataFrame

    Build the terminal node table from the loaded trees in `graphs`.

    Each column (`x1`, `x2`, …) is one terminal node. The number of terminals is taken
    from the `"V"` graph, and the preterminal edges of every tree are assumed to be in the
    same order, so that the k-th preterminal edge of each tree belongs to terminal k.
    The rows are blocks stacked on top of each other (`n_in` = number of inflow trees,
    `n_trees` = number of all trees):

    | Rows      | Content                                                              |
    |:----------|:---------------------------------------------------------------------|
    | `1+n_in`  | species ids: row 1 is the terminal `T_k`, then one row per inflow tree |
    | `n_trees` | coordinates `(x, y, z)` of the terminal, one row per tree            |
    | `n_trees` | tree id each coordinate row belongs to                               |
    | `1+n_in`  | flow ids: row 1 is `"Outflow"`, then one row per inflow tree         |
    | `1+n_in`  | flow values: row 1 is the sum of all outflow flows, then one row per inflow tree |
    | `1`       | terminal volume, `volume_geometry / n_terminals`                     |

    Inflow rows and coordinate rows follow the order of `tree_info.vascular_trees`.

    # Matching edges to node coordinates

    Only nodes have coordinates, so each preterminal edge takes the coordinates of one of
    its two nodes (the attribute convention is described in `ReadTreeFiles`). Edges in
    outflow trees were reversed earlier, so the leaf node is on a different side:

    - inflow trees: the target node, e.g. edge `(3, 6)` gets the coordinates of node 6;
    - outflow trees: the source node.

    `equality_edge_to_node` holds the index of that node in the edge tuple: `1` for the
    source, `2` for the target.
    """
    function prepare_tnodes_info(
        graphs::OrderedDict{String,DataFrame},
        tree_info,
    )::DataFrame
        @info "Preparation steps"
        # next lines are for being able to allocate vectors below
        # get the number of terminal nodes
        # not a good way, but for now its okay
        n_terminals::Int = length(
            selection_from_df(
                graphs["V"],
                (graphs["V"].ODE_groups .== groups.terminal, :ODE_groups),
            ),
        )
        # get number of inflows for this tree
        n_inflow::Int = length(tree_info.tree_components[:inflow_trees])
        # get number of all graphs (vascular trees)
        n_vascular_trees::Int = length(tree_info.tree_components[:inflow_trees]) + length(tree_info.tree_components[:outflow_trees])

        # === preallocating arrays which will store values/strings that are needed to write correct ODEs for terminal part
        # or to be able to understand to which element of the model values belongs
        flow_values::Array{Float64,2} = zeros(n_inflow + 1, n_terminals)
        flow_ids::Array{String,2} = fill("", n_inflow + 1, n_terminals)
        species_ids::Array{String,2} = fill("", n_inflow + 1, n_terminals)
        species_coordinates::Array{Tuple} = fill((0.0, 0.0, 0.0), n_vascular_trees, n_terminals)
        # to what tree part belong every species
        species_affiliations::Array{String,2} = fill("", n_vascular_trees, n_terminals)

        # where the inflow species, concentrations should be stored (row number in an array)
        # inflows start from the row two
        k_inflow::Int = 2

        p = Progress(length(values(graphs)); dt=0.5, color=:magenta)
        # iterate over all graphs' dataframes
        for (kg, graph) in enumerate(values(graphs))
            if graph[1, :is_inflow]
                # map edge with node by target id (see function's description)
                equality_edge_to_node = 2
                # get flow's values from this graph
                flow_values[k_inflow, :] .=
                    selection_from_df(graph, (graph.ODE_groups .== groups.preterminal, :flows))
                # get flow's ids from this graph
                flow_ids[k_inflow, :] .= selection_from_df(
                    graph,
                    (graph.ODE_groups .== groups.preterminal, :flow_ids),
                ) 
                # get species' ids from this graph
                species_ids[k_inflow, :] .= selection_from_df(
                    graph,
                    (graph.ODE_groups .== groups.preterminal, :species_ids),
                )
                # update inflows' iterator: next inflow values will go to the next row 
                k_inflow += 1
            else
                # map edge with node by source id (see function's description)
                equality_edge_to_node = 1
                # because of the ODE's struncture for the terminal part, we can sum outflows' flows,
                # but not store them separately for the place economy and for convinience
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
                        break
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
        tnodes_info_df = vcat(df_species_ids, df_species_coordinates, df_species_affiliations, df_flow_ids, df_flow_values)
        push!(tnodes_info_df, [volume_terminal for _ = 1:n_terminals])
        printl(tnodes_info_df)
        return tnodes_info_df
    end

end
