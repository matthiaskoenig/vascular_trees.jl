using VascularTrees
using VascularTrees.ProcessIndividualTrees: process_julia_graph
using VascularTrees.ProcessTerminalNodes: process_terminal_nodes
using VascularTrees.Simulation: run_simulations
using VascularTrees.Visualization: write_vtkhdf_files

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
    write_vtkhdf_files(tree_info, 1.0, 0.02)   # dt = 16.0 / 800
end

# in the julia terminal with activated environment run
# using PackageCompiler
# create_sysimage(
#     [
#         "VascularTrees",
#         ODE solving
#         "DifferentialEquations", "OrdinaryDiffEq", "Sundials", "SciMLBase", "SparseArrays",
#         tables / IO
#         "Arrow", "CSV", "DataFrames", "DataFramesMeta", "Tables", "Dictionaries", "HDF5",
#         utilities used inside the package
#         "Parameters", "Reexport", "ProgressMeter", "TimerOutputs",
#     ];
#     sysimage_path = "JuliaSysimage.so",
#     project = ".",
#     precompile_execution_file = "scripts/precompile_project.jl",
# )
