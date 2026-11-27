using VascularTrees
using VascularTrees.TreesProcessing:
    calculate_pressure_distribution

# Definitions of the available tree configurations (inflow / outflow trees)
const trees::TreeDefinitions = TreeDefinitions()

t_options = TreeOptions(
    n_nodes = [10],  #, 10, 100, 1000, 10000, 100000, 1000000, 10000000
    tree_configurations = [
        "Rectangle_quad",
        #"Rectangle_trio",
    ],
)
for tree_configuration ∈ t_options.tree_configurations
    for n_node ∈ t_options.n_nodes
        # fill the tree structure with the info for this tree
        tree_info =
            TreeStructure(; tree_configuration = tree_configuration, n_node = n_node, tree_components = trees.vascular_trees[tree_configuration])
        @info "Working on $(tree_info.graph_id)"
        (sol, unknowns) = calculate_pressure_distribution(tree_info; html_report = true)
        @show sol.retcode
    end
end
