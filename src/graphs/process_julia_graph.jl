module Process_Julia_Graph
    """
    Module for processing julia graph files from SyntheticVascularTrees.jl.
        Only works with "Rectangle_quad" and "Rectangle_trio". Otherwise "tree_definitions" must be modified in Utils.jl.

    Input: .lg and .csv files for every individual vessel tree (arterial, portal, etc.)

    Output: .arrow (table) for every individual vessel tree (arterial, portal, etc.) + 
        one table for terminal nodes which serve as connection points between the trees.

    Idea of this module: to get, to store, and to save all the information that we need for correct ODEs from the graph. 

    TODO: Optimize code
    """

    include("../utils.jl")
    import .Utils.Options: graph_options
    include("./processing_helpers.jl")

    include("./process_individual_trees.jl")
    import .Process_Individual_Trees: process_julia_graph
    include("./process_terminal_nodes.jl")
    import .Process_Terminal_Nodes: process_terminal_nodes

    # ============ Specify options
    g_options = graph_options(
        n_nodes = [10],  #750, 1000, 1250, 1500
        tree_ids = [
            "Rectangle_quad",
        ],
    )

    for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
        process_julia_graph(tree_id, n_node)
        process_terminal_nodes(tree_id, n_node)
    end

end
