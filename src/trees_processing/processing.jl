"""
Functions that convert the vessel trees written by the SyntheticVascularTrees generator into the
`.arrow` tables used by the simulations.

Data flow:
1. `process_julia_graph` (individual_trees.jl): for every individual tree (A, P, V, B, ...) the
   generator files are read (read_tree_files.jl), special edges are labeled and created
   (edge_labeling.jl), and `<GRAPH_DIR>/graphs/<tree>.arrow` is saved.
2. `process_terminal_nodes` (terminal_nodes.jl): the `.arrow` tables of all trees are loaded,
   the terminal nodes of the different trees are matched (closest_nodes.jl), and
   `<GRAPH_DIR>/graphs/T.arrow` is saved.

Step 1 must run before step 2.
"""
module TreesProcessing

    export process_julia_graph,
        process_terminal_nodes

    using DataFrames,
        Parameters,
        Arrow,
        ProgressMeter,
        CSV,
        NearestNeighbors,
        StaticArrays,
        Statistics,
        Hungarian,
        Combinatorics

    using ..Definitions:
        FlowDirections,
        ODE_groups

    using ..TableUtils:
        get_extended_vector,
        rows_to_tuples,
        selection_from_df,
        save_as_arrow,
        swap_parts,
        block_to_df,
        load_arrow

    # ODE group codes and inflow/outflow tree names, defined in definitions.jl
    const groups::ODE_groups = ODE_groups()
    const flow_direction = FlowDirections()

    # Volume of the liver cube in which the tree was grown (0.1 cm sides): [cm^3] -> [L].
    # It is split equally between all terminal nodes.
    const volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000

    include("read_tree_files.jl")
    include("trees/edge_labeling.jl")
    include("trees/individual_trees.jl")
    include("terminal_nodes/closest_nodes.jl")
    include("terminal_nodes/terminal_nodes.jl")

end
