"""
Write the simulation results as time-dependent VTKHDF files for ParaView.

Run: `julia --project scripts/create_vtk.jl` (writing starts as soon as the file is included).

What must be specified below:
1. t_options (tree configurations and numbers of nodes)
2. flow_scaling_factors and dt - must be the values of the simulations
   (dt = (tspan[2] - tspan[1]) / steps of the SimulationOptions in run_simulations.jl)

Inputs, in the graph directory
`results/julia_vessel_trees/<tree_configuration>/<tree_configuration>_<n_nodes>/`:
1. `graphs/<tree>.arrow` - written by build_arrow_graphs.jl
2. `simulations/<tree>[_<f>Q]_simulations_<dt>_dt.arrow` - written by run_simulations.jl

Outputs, in the `vtkhdf/` folder of the graph directory (created if needed):
1. `<tree>[_<f>Q]_<dt>_dt.vtkhdf` - one file per vascular tree (not for the terminal part T):
   lines of the tree with `Flow`, `Volume` and the time-dependent `Concentration` per edge.
   Open in ParaView ≥ 5.12; after rewriting a file that is open in ParaView, use
   File > Reload Files.
"""
module CreateVTK
    # === Imports ===
    using Revise
    using VascularTrees
    using VascularTrees.Visualization: write_vtkhdf_files

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

    # time step of the simulations [min], part of the simulation file names
    dt = 0.02


    # ============================== Run ==============================

    """
        create_all_vtks()

    Write the VTKHDF files for every combination of tree configuration, number of nodes and flow
    scaling factor from the options above. Order: configurations (outermost), number of nodes,
    flow scaling factors (innermost).
    """
    function create_all_vtks()
        # product varies its first iterator fastest, so the order is reversed here
        for (flow_scaling_factor, n_node, tree_configuration) ∈ Iterators.product(
            flow_scaling_factors,
            t_options.n_nodes,
            t_options.tree_configurations,
        )
            tree_info =
                TreeStructure(; tree_configuration = tree_configuration, n_node = n_node, tree_components = trees.vascular_trees[tree_configuration])
            @info "Working on $(tree_info.graph_id)"
            write_vtkhdf_files(tree_info, flow_scaling_factor, dt)
        end
    end

    create_all_vtks()

end
