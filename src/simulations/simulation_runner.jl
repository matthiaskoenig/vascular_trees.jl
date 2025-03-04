
module Simulation_runner
"""
Module that summarizes functions that load graph file in arrow format and run simulations.

Make sure that structure of vascular_trees.jl directory is as needed (somewhere it should be written)

What must/may be specified in code below:
1. must
    1.1. g_options (graph options)
    1.2. sim_options (simulation options)
2. may (if you whant to geet additional data and do not like default variants)
    2.1. benchmark options 
    2.2. solver options

Inputs:
1. Graph file (prepared information in arrow format (process_julia_graph.jl) from graph generated in SyntheticVascularTrees.jl)

Showed in the terminal:
1. Julia stuff
2. benchmark results if sim_options.benchmark=true

Outputs:
1. simulations.csv - if you do not do benchmarking and sim_options.save_simulations = true
2. running_times.csv - if sim_options.benchmark abd bench_options.save_running_times = true

TODO: run macros @timeit only when you need it - now it is so, but it is dirty
"""

include("../utils.jl")
import .Utils.Options: graph_options, simulations_options, benchmark_options, solver_options

include("simulation_helpers.jl")
import .Simulation_helpers: create_simulations

using InteractiveUtils

# ============ Specify options
g_options::graph_options = graph_options(
    n_nodes = [10],  #10
    tree_ids = [
        "Rectangle_quad",
        # "Rectangle_trio",
    ],
)

tspan = (0.0, 15.0)
sim_options::simulations_options = simulations_options(
    tspan = tspan,
    tpoints = range(tspan[1], stop = tspan[2], length = 30),
    save_simulations = false,
    benchmark = false,
)

sol_options::solver_options = solver_options()

if sim_options.benchmark

    import .Utils.Options: benchmark_options

    # do not write anything here in brackets if you are okay with default variant
    bench_options::benchmark_options = benchmark_options(save_running_times=false)
    create_simulations(
        g_options,
        sim_options,
        bench_options,
        sol_options,
    )

else
    create_simulations(
        g_options,
        sim_options,
        sol_options,
    )
end

end
