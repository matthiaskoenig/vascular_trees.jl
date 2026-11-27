using DataFrames
using VascularTrees.TreesProcessing: process_julia_graph, process_terminal_nodes
using VascularTrees.TableUtils: load_arrow

# End-to-end run on the smallest tree, in a temporary copy of the input files,
# so that results/ is not touched.
GRAPH_DIR = joinpath(mktempdir(), "Rectangle_quad_10")
cp(joinpath(DATA_DIR, "Rectangle_quad_10"), GRAPH_DIR)

tree_info = TreeStructure(;
    tree_configuration = "Rectangle_quad",
    n_node = 10,
    tree_components = TreeDefinitions().vascular_trees["Rectangle_quad"],
    GRAPH_DIR = GRAPH_DIR,
)
m_options = MatchingOptions(
    k_greedy = 3,
    max_iter = 10,
    tol = 1e-8,
    algorithm = AuctionAssignment(),
    compute_bound = false,
)
flow_scaling_factor = 1.0
find_closest_tnodes = false

process_julia_graph(tree_info)
process_terminal_nodes(tree_info, find_closest_tnodes, m_options)
for tree ∈ [tree_info.vascular_trees; "T"]
    @test isfile(joinpath(GRAPH_DIR, "graphs", "$(tree).arrow"))
end