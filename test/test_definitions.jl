trees = TreeDefinitions()

@test Set(keys(trees.vascular_trees)) == Set(["Rectangle_trio", "Rectangle_quad"])

tree_info = TreeStructure(;
    tree_configuration = "Rectangle_quad",
    n_node = 10,
    tree_components = trees.vascular_trees["Rectangle_quad"],
)
@test tree_info.graph_id == "Rectangle_quad_10"
@test sort(tree_info.vascular_trees) == ["A", "B", "P", "V"]
@test endswith(tree_info.GRAPH_DIR, joinpath("Rectangle_quad", "Rectangle_quad_10"))

# GRAPH_DIR can point outside of the results folder (used by test_pipeline.jl)
custom_dir = mktempdir()
tree_info = TreeStructure(;
    tree_configuration = "Rectangle_trio",
    n_node = 10,
    tree_components = trees.vascular_trees["Rectangle_trio"],
    GRAPH_DIR = custom_dir,
)
@test tree_info.GRAPH_DIR == custom_dir
@test sort(tree_info.vascular_trees) == ["A", "P", "V"]

# dt is derived from tspan and steps
sim_options = SimulationOptions(tspan = (0.0, 16.0), steps = 800, save_simulations = true, benchmark = false)
@test sim_options.dt ≈ 0.02
