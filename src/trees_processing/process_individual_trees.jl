"""
    ProcessIndividualTrees

Converts every individual vessel tree (e.g. arterial `A`, portal `P`, venous `V`, biliary `B`)
from the graph files written by the tree generator into one `.arrow` table that contains
everything the ODE model and the terminal-nodes processing need.

# Input (per tree, in `<GRAPH_DIR>/julia/`)
- `<tree>.lg`        - list of edges (graph skeleton).
- `<tree>_edges.csv` - edge attributes (radius, flow, length, ...), stored per target node.
- `<tree>_nodes.csv` - node attributes (id, x/y/z coordinates).

# Output
- `<GRAPH_DIR>/graphs/<tree>.arrow` - one table per tree, see [`build_arrow_table`](@ref).

# Main steps
1. read the files into an edges dataframe and a nodes dataframe;
2. add ids, edge labels, terminal and marginal edges, ODE groups;
3. reverse edges of outflow trees;
4. find predecessors and successors of every edge;
5. collect everything in one table and save it.

TODO: ids of outflow edges are wrong: they are built before the edges are reversed, so they are
`source_target` of the original direction (`target_source` of the final one).
"""
module ProcessIndividualTrees

    export process_julia_graph

    using ..Definitions: flow_directions, ODE_groups
    using ..ReadTreeFiles: read_edges, read_nodes
    using ..TableUtils: get_extended_vector, create_tuples_from_dfrows, selection_from_df, save_as_arrow

    using ..ProcessingHelpers: label_special_edges!, create_special_edges!

    using DataFrames, InteractiveUtils

    # ODE group codes and inflow/outflow tree names, defined in definitions.jl
    const groups::ODE_groups = ODE_groups()
    const flow_direction = flow_directions()

    """
        process_julia_graph(tree_info)

    Process every vessel tree listed in `tree_info.vascular_trees` and save one `.arrow` file per tree.

    # Arguments
    - `tree_info`: a `Tree_structure` (provides `vascular_trees` and `GRAPH_DIR`).
    """
    function process_julia_graph(tree_info)
        @info "Processing individual trees..."
        for vascular_tree ∈ tree_info.vascular_trees
            process_individual_tree(tree_info.GRAPH_DIR, vascular_tree)
        end
    end

    """
        process_individual_tree(GRAPH_DIR, vascular_tree)

    Full pipeline for one vessel tree: read its files, add the information needed for the ODE
    system, build the output table and save it to `<GRAPH_DIR>/graphs/<vascular_tree>.arrow`.
    """
    function process_individual_tree(GRAPH_DIR::String, vascular_tree::String)
        GRAPH_PATH, EDGES_PATH, NODES_PATH = paths_initialization(GRAPH_DIR, vascular_tree)
        # graph skeleton and edge attributes are joined in one dataframe
        edges_df, nodes_df = read_graph(GRAPH_PATH, EDGES_PATH, NODES_PATH)
        add_graph_characteristics!(edges_df, vascular_tree)
        graph = build_arrow_table(edges_df, nodes_df, vascular_tree)
        save_as_arrow(graph, joinpath(GRAPH_DIR, "graphs/$(vascular_tree).arrow"))
    end

    #=================================================================================================================================#
    """
        paths_initialization(GRAPH_DIR, vascular_tree) -> (GRAPH_PATH, EDGES_PATH, NODES_PATH)

    Return paths to the `.lg`, `_edges.csv` and `_nodes.csv` files of `vascular_tree`
    (all in `<GRAPH_DIR>/julia/`).
    """
    function paths_initialization(
        GRAPH_DIR::String,
        vascular_tree::String,
    )::Tuple{String,String,String}
        GRAPH_PATH::String = joinpath(GRAPH_DIR, "julia/$(vascular_tree).lg")
        EDGES_PATH::String = joinpath(GRAPH_DIR, "julia/$(vascular_tree)_edges.csv")
        NODES_PATH::String = joinpath(GRAPH_DIR, "julia/$(vascular_tree)_nodes.csv")

        return GRAPH_PATH, EDGES_PATH, NODES_PATH
    end

    """
        read_graph(GRAPH_PATH, EDGES_PATH, NODES_PATH) -> (edges_df, nodes_df)

    Read the graph files of one tree.

    # Returns
    - `edges_df`: one row per edge; `source_id`, `target_id` joined with the edge
      attributes (`leaf`, `radius`, `flows` [L/min], `length`, `pressure_drop`, `volumes` [L]).
    - `nodes_df`: one row per node; `ids`, `x`, `y`, `z`.

    See `read_edges` / `read_nodes` in helpers.jl for renaming and unit conversion.
    """
    function read_graph(
        GRAPH_PATH::String,
        EDGES_PATH::String,
        NODES_PATH::String,
    )::Tuple{DataFrame,DataFrame}
        edges_df::DataFrame = read_edges(GRAPH_PATH, EDGES_PATH)
        nodes_df::DataFrame = read_nodes(NODES_PATH)

        return edges_df, nodes_df
    end

    """
        add_graph_characteristics!(edges_df, vascular_tree)

    Add to the edges dataframe everything the ODE model needs, in this order:
    1. `species_ids`, `flow_ids`, `volume_ids` - names `<tree>_<source>_<target>`, `Q_<source>_<target>`,
       `V_<source>_<target>`, used as column names of simulation results;
    2. boolean labels `preterminal`, `start`, `terminal`;
    3. new rows: terminal self-edges `(n, n)` and the marginal edge `(0, start_node)`;
    4. `ODE_group` of every edge;
    5. `index` - row position of every edge.

    The order matters: terminal edges can only be created after preterminal edges are labelled,
    and ODE groups / indices must be assigned after all rows exist.
    """
    function add_graph_characteristics!(edges_df::DataFrame, vascular_tree::String)
        # ids are used to name the columns of simulation results and to map every
        # flow/volume value back to its edge
        transform!(
            edges_df,
            [:source_id, :target_id] =>
                ByRow(
                    (source_id, target_id) -> (
                        "$(vascular_tree)_$(source_id)_$(target_id)",
                        "Q_$(source_id)_$(target_id)",
                        "V_$(source_id)_$(target_id)",
                    ),
                ) => [:species_ids, :flow_ids, :volume_ids],
        )
        # terminal edges don't exist yet, so `terminal` is false for every edge here
        label_special_edges!(edges_df)
        # one self-edge per terminal node + the marginal edge, where the dose is applied
        create_special_edges!(edges_df)
        # one integer per edge lets the ODE function choose the right equation quickly
        assign_ODE_group!(edges_df)
    end

    """
        build_arrow_table(edges_df, nodes_df, vascular_tree) -> DataFrame

    Build the table that is saved as `<vascular_tree>.arrow`.

    The table has one row per edge (including terminal and marginal edges). Arrow can only store
    columns of equal length, so columns that are shorter are padded with `missing`.

    # Columns - one value per edge
    - `all_edges`: `(source_id, target_id)`; outflow trees are already reversed.
    - `flows`, `volumes`: flow [L/min] and volume [L] of the edge.
    - `species_ids`, `flow_ids`, `volume_ids`: names of the edge's variables.
    - `ODE_groups`: see `ODE_groups` in definitions.jl.
    - `pre_elements`, `post_elements`: row indices of predecessor / successor edges.

    # Columns - padded with `missing`
    - `vascular_tree_id`, `is_inflow`: one value (first row).
    - `nodes_ids`, `nodes_coordinates`: one value per node, coordinates as `(x, y, z)`.
    - `terminal_edges`, `start_edge`, `preterminal_edges`: `(source_id, target_id)` of those edges.

    Note: modifies `edges_df` (edges of outflow trees are reversed, see [`prepare_graph_info`](@ref)).
    """
    function build_arrow_table(
        edges_df::DataFrame,
        nodes_df::DataFrame,
        vascular_tree::String
    )::DataFrame
        graph_info::NamedTuple =
            prepare_graph_info(edges_df, nodes_df, vascular_tree)
        df_length::Int = length(graph_info.all_edges)

        graph = DataFrame(
            vascular_tree_id = get_extended_vector(vascular_tree, df_length),
            is_inflow = get_extended_vector(graph_info.is_inflow, df_length),
            nodes_ids = get_extended_vector(nodes_df.ids, df_length),
            nodes_coordinates = get_extended_vector(graph_info.nodes_coordinates, df_length),
            all_edges = graph_info.all_edges,
            terminal_edges = get_extended_vector(graph_info.terminal_edges, df_length),
            start_edge = get_extended_vector(graph_info.start_edge, df_length),
            preterminal_edges = get_extended_vector(graph_info.preterminal_edges, df_length),
            flows = edges_df.flows,
            volumes = edges_df.volumes,
            species_ids = edges_df.species_ids,
            flow_ids = edges_df.flow_ids,
            volume_ids = edges_df.volume_ids,
            ODE_groups = graph_info.ODE_groups,
            pre_elements = graph_info.pre_elements,
            post_elements = graph_info.post_elements,
        )

        return graph
    end

    #=================================================================================================================================#
    """
        assign_ODE_group!(edges_df)

    Add column `ODE_group`. Checked in this order (first match wins):
    preterminal → `groups.preterminal`, terminal → `groups.terminal`,
    marginal (`source_id == 0`) → `groups.marginal`, everything else → `groups.other`.
    """
    function assign_ODE_group!(edges_df::DataFrame)
        """Function that adds to the dataframe with edges column which indicate ODE group for each edge"""
        edges_df.ODE_group .=
        ifelse.(
            edges_df[:, :preterminal] .== true,
            groups.preterminal,
            ifelse.(
                edges_df[:, :terminal] .== true,
                groups.terminal,
                ifelse.(
                    edges_df[:, :source_id] .== 0,
                    groups.marginal,
                    groups.other,
                ),
            ),
        )
    end


    """
        prepare_graph_info(edges_df, nodes_df, vascular_tree) -> NamedTuple

    Extract from both dataframes the vectors needed by [`build_arrow_table`](@ref).

    For outflow trees the edges are reversed **in place** (`source_id` and `target_id` columns
    are swapped), so that every tree is described in the direction of blood flow.

    # Returns
    NamedTuple with `is_inflow`, `nodes_coordinates`, `all_edges`, `terminal_edges`, `start_edge`,
    `preterminal_edges`, `ODE_groups`, `pre_elements`, `post_elements`.
    Edges and coordinates are tuples: `(source_id, target_id)` / `(x, y, z)`.
    """
    function prepare_graph_info(
        edges_df::DataFrame,
        nodes_df::DataFrame,
        vascular_tree::String,
    )::NamedTuple
        is_inflow::Bool = in(vascular_tree, flow_direction.inflow_trees)
        # outflow trees: swap source and target so that edges follow the blood flow
        (!is_inflow) &&
            (rename!(edges_df, [:source_id => :target_id, :target_id => :source_id]))

        edges = selection_from_df(edges_df, (:, [:source_id, :target_id]))
        terminals = selection_from_df(
            edges_df,
            (edges_df.terminal .== true, [:source_id, :target_id]),
        )
        start = selection_from_df(
            edges_df,
            (edges_df.start .== true, [:source_id, :target_id]),
        )
        preterminals = selection_from_df(
            edges_df,
            (edges_df.preterminal .== true, [:source_id, :target_id]),
        )
        ODE_groups = selection_from_df(edges_df, (:, :ODE_group))
        nodes_coord = selection_from_df(nodes_df, (:, [:x, :y, :z]))

        # tuples are what is stored in the .arrow file
        all_edges, terminal_edges, start_edge, preterminal_edges, nodes_coordinates =
            map(create_tuples_from_dfrows, (edges, terminals, start, preterminals, nodes_coord))

        pre_elements, post_elements = get_pre_postelements(all_edges)

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

    """
        get_pre_postelements(edges) -> (pre_elements, post_elements)

    For every edge `(s, t)` find:
    - predecessors: edges that end in `s`;
    - successors: edges that start in `t`.
    The edge itself is excluded (terminal self-edges `(n, n)` would otherwise match themselves).

    They tell the ODE function where each edge gets its inflow from and
    where it sends its outflow.
    """
    function get_pre_postelements(edges::Vector{Tuple{Int,Int}})
        incoming = Dict{Int,Vector{Int}}()   # node id => indices of edges ending at it
        outgoing = Dict{Int,Vector{Int}}()   # node id => indices of edges starting at it
        for (i, (s, t)) in enumerate(edges)
            push!(get!(outgoing, s, Int[]), i)
            push!(get!(incoming, t, Int[]), i)
        end
        # filter(!=(i)) excludes the edge itself (terminal edges have source == target)
        pre  = [filter(!=(i), get(incoming, s, Int[])) for (i, (s, _)) in enumerate(edges)]
        post = [filter(!=(i), get(outgoing, t, Int[])) for (i, (_, t)) in enumerate(edges)]
        return pre, post
    end

end
