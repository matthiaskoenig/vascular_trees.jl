using DataFrames
using VascularTrees.TreesProcessing: process_julia_graph, process_terminal_nodes
using VascularTrees.Simulation: run_simulations
using VascularTrees.Visualization: write_vtkhdf_files
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
sim_options = SimulationOptions(
    tspan = (0.0, 1.0),
    steps = 10,
    save_simulations = true,
    benchmark = false,
    show_progress = false,
)
flow_scaling_factor = 1.0
find_closest_tnodes = false

process_julia_graph(tree_info)
process_terminal_nodes(tree_info, find_closest_tnodes, m_options)
for tree ∈ [tree_info.vascular_trees; "T"]
    @test isfile(joinpath(GRAPH_DIR, "graphs", "$(tree).arrow"))
end

run_simulations(
    tree_info,
    sim_options,
    SolverOptions(),
    (dense = false, save_everystep = false, progress = false),
    flow_scaling_factor,
    BenchmarkOptions(),
)
for tree ∈ [tree_info.vascular_trees; "T"]
    sim_path = joinpath(GRAPH_DIR, "simulations", "$(tree)_simulations_$(sim_options.dt)_dt.arrow")
    @test isfile(sim_path)
    concentrations = Matrix(load_arrow(sim_path)[:, 2:end])
    @test size(concentrations, 2) == sim_options.steps + 1
    @test all(isfinite, concentrations)
    @test all(>=(-1e-10), concentrations)
end

write_vtkhdf_files(tree_info, flow_scaling_factor, sim_options.dt)
for tree ∈ tree_info.vascular_trees
    @test isfile(joinpath(GRAPH_DIR, "vtkhdf", "$(tree)_$(sim_options.dt)_dt.vtkhdf"))
end
