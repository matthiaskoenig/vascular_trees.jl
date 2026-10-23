using VascularTrees
using VascularTrees.ProcessIndividualTrees: process_julia_graph
using VascularTrees.ProcessTerminalNodes: process_terminal_nodes
using VascularTrees.Simulation: run_simulations

trees = TreeDefinitions()
for cfg in ["Rectangle_quad"]           
    tree_info = TreeStructure(; tree_configuration = cfg, n_node = 10,
                              tree_components = trees.vascular_trees[cfg])
    process_julia_graph(tree_info)
    process_terminal_nodes(tree_info)
    run_simulations(
        tree_info,
        SimulationOptions(tspan = (0.0, 16.0), 
        steps = 800,
        save_simulations = true, 
        benchmark = false
        ),
        SolverOptions(),
        (dense = false, save_everystep = false, progress = false),
        1.0,
        BenchmarkOptions()
    )
end

# in the julia terminal with activated environment run
# using PackageCompiler
# create_sysimage(
#     [
#         "VascularTrees",
#         # ODE solving
#         "DifferentialEquations", "OrdinaryDiffEq", "Sundials", "SciMLBase", "SparseArrays",
#         # tables / IO
#         "Arrow", "CSV", "DataFrames", "DataFramesMeta", "Tables", "Dictionaries",
#         # utilities used inside the package
#         "Parameters", "Reexport", "ProgressMeter", "TimerOutputs",
#     ];
#     sysimage_path = "VascularTrees_sysimage.so",
#     project = ".",
#     precompile_execution_file = "scripts/precompile_project.jl",
# )
