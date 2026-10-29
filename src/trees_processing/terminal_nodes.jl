"""
Module for collecting information about terminal nodes from the files of all individual
vessel trees (arterial, portal, etc.).

Input: .arrow (table) for every individual vessel tree (arterial, portal, etc.).

Output: `T.arrow` (table) with information about terminal nodes.

Idea of this module: to collect and save the information about terminal nodes that is needed
    for the terminal part of the ODEs. These nodes are the connection points between the
    individual vessel trees.

Workflow:
1. [`process_terminal_nodes`](@ref) - entry point: load, process, save.
2. [`collect_graphs_information!`](@ref) - load the `.arrow` table of every tree.
3. [`prepare_tnodes_info`](@ref) - fill the blocks of the output table, tree by tree.
    - [`get_terminal_coordinates`](@ref) - coordinates of every terminal in one tree.
    - [`assemble_tnodes_table`](@ref) - stack the blocks into one `DataFrame`.

TODO: Flow affiliations are wrong.

TODO: prepare_tnodes_info assumes that every tree lists its preterminal edges in the same order,
    so that the k-th preterminal edge of each tree belongs to terminal k (column k of the output).
    This is not checked. If it does not hold, values of different terminals get mixed up
    (possible cause of the wrong flow affiliations).
"""
module ProcessTerminalNodes

    export process_terminal_nodes

    using ..Definitions: ODE_groups
    using ..TableUtils: 
        selection_from_df, 
        save_as_arrow, 
        get_extended_vector,
        block_to_df,
        load_arrow

    using DataFrames,
        Arrow,
        ProgressMeter

    using CSV

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

        # tree name => its table; the order of the output rows comes from
        # tree_info.vascular_trees, not from this (unordered) Dict
        graphs = Dict{String,DataFrame}()
        # load the tables of all trees (graphs) that belong to this tree configuration
        collect_graphs_information!(graphs, tree_info)

        tnodes_info_df = prepare_tnodes_info(graphs, tree_info)
        save_as_arrow(tnodes_info_df, joinpath(tree_info.GRAPH_DIR, "graphs/T.arrow"))
        # CSV.write(joinpath(tree_info.GRAPH_DIR, "graphs/T.csv"),tnodes_info_df)
    end

    """
        collect_graphs_information!(graphs::Dict{String,DataFrame}, tree_info)

    Load the `.arrow` table of every tree in `tree_info.vascular_trees` into `graphs`,
    keyed by the tree name (e.g. `"A"`, `"V"`). Only the columns needed for the terminal
    nodes are loaded.
    """
    function collect_graphs_information!(graphs::Dict{String,DataFrame}, tree_info)
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
        prepare_tnodes_info(graphs::Dict{String,DataFrame}, tree_info) -> DataFrame

    Build the terminal node table from the loaded trees in `graphs`.

    The output has one row per terminal node. The number of terminals is taken from the `"V"`
    graph, and the preterminal edges of every tree are assumed to be in the same order, so that
    the k-th preterminal edge of each tree belongs to terminal k (row k).

    Inside this function the data are collected in blocks of size `rows × n_terminals`
    (`n_in` = number of inflow trees, `n_trees` = number of all trees), column k of a block
    belongs to terminal k. [`assemble_tnodes_table`](@ref) turns each block into columns of the
    output, see there for the column names.

    | Block                  | Rows      | Content                                                              |
    |:-----------------------|:----------|:---------------------------------------------------------------------|
    | `species_ids`          | `1+n_in`  | row 1 is the terminal `T_k`, then one row per inflow tree            |
    | `species_coordinates`  | `n_trees` | coordinates `(x, y, z)` of the terminal, one row per tree            |
    | `species_affiliations` | `n_trees` | tree id each coordinate row belongs to                               |
    | `flow_ids`             | `1+n_in`  | row 1 is `"Outflow"`, then one row per inflow tree                   |
    | `flow_values`          | `1+n_in`  | row 1 is the sum of all outflow flows, then one row per inflow tree  |

    Inflow rows and coordinate rows follow the order of `tree_info.vascular_trees`.
    Coordinates come from [`get_terminal_coordinates`](@ref). A terminal whose leaf node
    is not found keeps the coordinates `(0.0, 0.0, 0.0)` and the tree id `""`.
    """
    function prepare_tnodes_info(
        graphs::Dict{String,DataFrame},
        tree_info,
    )::DataFrame
        @info "Preparation steps"
        # === Table sizes, needed to preallocate the arrays below
        # Number of terminal nodes = number of columns of the output table.
        # TODO: taken from the "V" tree only; assumes every configuration has a "V" tree.
        n_terminals::Int = count(==(groups.terminal), graphs["V"].ODE_groups)
        n_inflow::Int = length(tree_info.tree_components[:inflow_trees])
        n_vascular_trees::Int = length(tree_info.vascular_trees)

        # === Preallocated blocks of the output table (row layout: see the docstring)
        # Column k belongs to terminal k in every block.
        # Blocks with `n_inflow + 1` rows: row 1 = the terminal itself / all outflow trees,
        # rows 2, 3, ... = one row per inflow tree.
        flow_values::Array{Float64,2} = zeros(n_inflow + 1, n_terminals)
        flow_ids::Array{String,2} = fill("", n_inflow + 1, n_terminals)
        flow_ids[1, :] .= "Outflow"
        species_ids::Array{String,2} = fill("", n_inflow + 1, n_terminals)
        species_ids[1, :] .= ["T_$n_terminal" for n_terminal = 1:n_terminals]
        # Blocks with one row per tree: coordinates of the terminal in that tree,
        # and the id of the tree the coordinates come from.
        species_coordinates::Array{NTuple{3,Float64}} = fill((0.0, 0.0, 0.0), n_vascular_trees, n_terminals)
        species_affiliations::Array{String,2} = fill("", n_vascular_trees, n_terminals)

        # Row that the next inflow tree writes to (row 1 is reserved for the outflow trees).
        k_inflow::Int = 2

        p = Progress(length(values(graphs)); dt=0.5, color=:magenta)
        # The order of tree_info.vascular_trees fixes the order of the rows.
        for (tree_row, tree_id) in enumerate(tree_info.vascular_trees)
            graph = graphs[tree_id]
            # rows of this tree's table that are preterminal edges (one per terminal)
            preterminal = graph.ODE_groups .== groups.preterminal
            if graph[1, :is_inflow]
                # the leaf of an inflow edge is its target node (see get_terminal_coordinates)
                leaf_node_index = 2
                # each inflow tree gets its own row with flows, flow ids and species ids
                flow_values[k_inflow, :] .=
                    selection_from_df(graph, (preterminal, :flows))
                flow_ids[k_inflow, :] .= selection_from_df(
                    graph,
                    (preterminal, :flow_ids),
                )
                species_ids[k_inflow, :] .= selection_from_df(
                    graph,
                    (preterminal, :species_ids),
                )
                k_inflow += 1
            else
                # outflow edges were reversed, so their leaf is the source node
                leaf_node_index = 1
                # The terminal ODEs only need the total outflow of each terminal, so the
                # flows of all outflow trees are added together in row 1.
                flow_values[1, :] .+=
                    selection_from_df(graph, (preterminal, :flows))
            end

            tree_name = tree_id #graph[1, :vascular_tree_id]

            # coordinates of every terminal in this tree; `found` marks the terminals
            # whose leaf node was found, only those get the tree id
            species_coordinates[tree_row, :], found = get_terminal_coordinates(
                graph.nodes_ids,
                graph.nodes_coordinates,
                graph.all_edges[preterminal],
                leaf_node_index
            )
            species_affiliations[tree_row, found] .= tree_name

            next!(p)
        end

        tnodes_info_df = assemble_tnodes_table(
            species_ids, 
            species_coordinates, 
            species_affiliations, 
            flow_ids, 
            flow_values, 
            n_terminals
        )
        return tnodes_info_df
    end

    """
        get_terminal_coordinates(nodes_ids, nodes_coordinates, preterminal_edges, leaf_node_index)
            -> (leaf_coordinates, found)

    Coordinates of every terminal in one tree: element k belongs to the k-th preterminal
    edge, i.e. to terminal k.

    Only nodes have coordinates, so each preterminal edge takes the coordinates of one of
    its two nodes, its leaf node (the attribute convention is described in `ReadTreeFiles`).
    Edges in outflow trees were reversed earlier, so the leaf node is on a different side:

    - inflow trees: the target node, e.g. edge `(3, 6)` gets the coordinates of node 6;
    - outflow trees: the source node.

    `leaf_node_index` is the index of the leaf node in the edge tuple: `1` for the source,
    `2` for the target.

    Takes the table columns, not the whole `DataFrame`, so that Julia compiles this
    function for the actual column types (function barrier).

    Returns the coordinates and a `BitVector` `found`. If the leaf node of an edge is not
    found, its coordinates stay `(0.0, 0.0, 0.0)` and `found` is `false`.
    """
    function get_terminal_coordinates(
        nodes_ids,
        nodes_coordinates,
        preterminal_edges,
        leaf_node_index,
    )
        leaf_coordinates::Vector{NTuple{3,Float64}} = fill((0.0, 0.0, 0.0), length(preterminal_edges))
        found = falses(length(preterminal_edges))
        # Node id => (x, y, z). Built once per tree so that each edge finds its node
        # directly instead of searching all nodes. Node columns are padded with
        # `missing` at the end; those rows are skipped.
        coordinates_by_node = Dict(
            id => coordinates
            for (id, coordinates) in zip(nodes_ids, nodes_coordinates)
            if !ismissing(id)
        )

        # k-th preterminal edge -> coordinates of terminal k as seen from this tree
        for (kg, edge) in enumerate(preterminal_edges)
            coordinates = get(coordinates_by_node, edge[leaf_node_index], nothing)
            # leaf node not found: keep the default (0, 0, 0), found[kg] stays false
            isnothing(coordinates) && continue
            leaf_coordinates[kg] = coordinates
            found[kg] = true
        end
        
        return leaf_coordinates, found
    end

    """
        assemble_tnodes_table(species_ids, species_coordinates, species_affiliations,
            flow_ids, flow_values, n_terminals) -> DataFrame

    Convert the blocks filled by [`prepare_tnodes_info`](@ref) into the output table:
    one row per terminal, every column has one concrete type (no `Any`).

    # Columns
    | Columns                        | Type                | Content                                               |
    |:-------------------------------|:--------------------|:------------------------------------------------------|
    | `species_id_1`                 | `String`            | terminal id `T_k`                                     |
    | `species_id_2`, …              | `String`            | species of the terminal in each inflow tree           |
    | `coordinates_1`, …             | `NTuple{3,Float64}` | `(x, y, z)` of the terminal in each tree              |
    | `affiliation_1`, …             | `String`            | tree id of each `coordinates_*` (`""` if not found)   |
    | `flow_id_1`                    | `String`            | `"Outflow"`                                           |
    | `flow_id_2`, …                 | `String`            | flow id of the terminal in each inflow tree           |
    | `flow_value_1`                 | `Float64`           | sum of the flows of all outflow trees                 |
    | `flow_value_2`, …              | `Float64`           | flow of the terminal in each inflow tree              |
    | `volume`                       | `Float64`           | terminal volume, `volume_geometry / n_terminals`      |

    The index in `species_id_*`, `flow_id_*` and `flow_value_*` runs over inflow trees (after the
    first column), in `coordinates_*` and `affiliation_*` over all trees, in the order of
    `tree_info.vascular_trees`.

    To get a block back in its `rows * n_terminals` form:
    `permutedims(Matrix(df[:, r"^coordinates_"]))`.
    """
    function assemble_tnodes_table(species_ids, species_coordinates, species_affiliations,
                                flow_ids, flow_values, n_terminals)
        return hcat(
            block_to_df(species_ids, "species_id"),          # species_id_1 = T_k, then inflow species
            block_to_df(species_coordinates, "coordinates"), # one column per tree
            block_to_df(species_affiliations, "affiliation"),
            block_to_df(flow_ids, "flow_id"),                # flow_id_1 = "Outflow"
            block_to_df(flow_values, "flow_value"),          # flow_value_1 = summed outflow
            DataFrame(volume = fill(volume_geometry / n_terminals, n_terminals)),
        )
    end

end
