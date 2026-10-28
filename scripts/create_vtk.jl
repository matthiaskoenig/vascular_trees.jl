module CreateVTK
    # === Imports ===
    using Revise
    using VascularTrees
    using VascularTrees.Visualization: write_vtkhdf

    # === Logging ===
    # colored log prefixes; debug = true also shows the per-subsystem stages of VascularTrees
    setup_logging(debug = false)

    # Definitions of the available tree configurations (inflow / outflow trees)
    const trees::TreeDefinitions = TreeDefinitions()


    # ============================ Options ============================

    # === Graph options ===
    # options for graph, i.e., number of nodes and type of tree
    t_options = TreeOptions(
        n_nodes = [10],  #, 100, 1000, 10000, 100000, 1000000
        tree_configurations = [
            "Rectangle_quad",
            #"Rectangle_trio",
        ],
    )

    # factors by which the flow is scaled, 1.0 is the original flow
    flow_scaling_factors = [1.0] # 1.0/16, 1.0/8, 1.0/4, 1.0/2, 1.0, 1.0*2, 1.0*4, 1.0*8, 1.0*16

    dt = 0.02


    function create_vtks()
        # product varies its first iterator fastest, so the order is reversed here
        for (flow_scaling_factor, n_node, tree_configuration) ∈ Iterators.product(
            flow_scaling_factors,
            t_options.n_nodes,
            t_options.tree_configurations,
        )
            tree_info =
                TreeStructure(; tree_configuration = tree_configuration, n_node = n_node, tree_components = trees.vascular_trees[tree_configuration])
            @info "Working on $(tree_info.graph_id)"

            for graph_subsystem in tree_info.vascular_trees
                print(graph_subsystem)
                ARROW_PATH::String = normpath(joinpath(tree_info.GRAPH_DIR, "graphs/$(graph_subsystem).arrow"))
                tree_df = load_arrow(ARROW_PATH)

                SIMULATIONS_PATH::String = normpath(joinpath(tree_info.GRAPH_DIR, "simulations/$(graph_subsystem)_simulations_$(dt)_dt.arrow")) 
                simulations_df = load_arrow(SIMULATIONS_PATH)

                VTK_DIR = normpath(joinpath(tree_info.GRAPH_DIR, "vtkhdf"))
                mkpath(VTK_DIR)
                VTK_PATH = normpath(joinpath(VTK_DIR, "$(graph_subsystem).vtkhdf"))

                write_vtkhdf(
                    tree_df,
                    simulations_df,
                    VTK_PATH
                )
            end
        end
    end

    create_vtks() 

end