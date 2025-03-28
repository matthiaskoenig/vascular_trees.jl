
module Simulation_Runner
"""
Module that summarizes functions that load graph file in arrow format and run simulations.

Make sure that structure of vascular_trees.jl directory is as needed (somewhere it should be written)

What must/may be specified in code below:
1. must
    1.1. g_options (graph options)
    1.2. sim_options (simulation options)
2. may (if you want to geet additional data and do not like default variants)
    2.1. benchmark options 
    2.2. solver options
    2.3. additional solver options (these arguments are not mandatory)

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
import .Utils: JULIA_RESULTS_DIR, MODEL_PATH 
import .Utils.Definitions: tree_definitions
import .Utils.Benchmarking: save_times_as_csv
include("../models/julia_from_jgraph.jl")
include("simulation_helpers.jl")
import .Simulation_Helpers: run_simulations


# Already specified in utils.jl
const trees::tree_definitions = tree_definitions()

# using InteractiveUtils

# === Graph options ===
# options for graph, i.e., number of nodes and type of tree
g_options = graph_options(
    n_nodes = [10, 100, 1000, 10000, 100000],  #10
    tree_configurations = [
        "Rectangle_quad",
        # "Rectangle_trio",
    ],
)

flow_scaling_factors = [1.0] # 1.0/16, 1.0/8, 1.0/4, 1.0/2, 1.0, 1.0*2, 1.0*4, 1.0*8, 1.0*16

# === Simulation options ===
sim_options = simulations_options(
    tspan = (0.0, 16.0),  # [min]
    steps = 800.0,
    save_simulations = true,
    benchmark = false,
)

# === ODE Solver options ===
# integrator, tolerances
sol_options = solver_options()
# additional integrator arguments
# these arguments are not mandatory
additional_sol_options::NamedTuple =
    (dense = false, save_everystep = false, progress = true)

# === Benchmark options ===
# do not write anything here in brackets if you are okay with default variant
import .Utils.Options: benchmark_options
bench_options = benchmark_options(save_running_times = false)

# Basic information about the tree that differs between its types (Rectangle_quad, trio, etc.)
# and which is used repeatedly in simulations
# DO NOT CHANGE
Base.@kwdef struct Tree_structure
    tree_configuration::String
    n_node::Integer
    graph_id::String = "$(tree_configuration)_$(n_node)"
    tree_components::Dict{Symbol,Vector{String}} = trees.vascular_trees[tree_configuration]
    vascular_trees::Vector{String} = reduce(vcat, values(tree_components))
    GRAPH_DIR::String = normpath(
        joinpath(@__FILE__, "../../..", JULIA_RESULTS_DIR, tree_configuration, graph_id),
    )
end

for tree_configuration ∈ g_options.tree_configurations
    for n_node ∈ g_options.n_nodes
        for flow_scaling_factor in flow_scaling_factors
            tree_info =
                Tree_structure(; tree_configuration = tree_configuration, n_node = n_node)
            run_simulations(
                tree_info,
                sim_options,
                sol_options,
                additional_sol_options,
                flow_scaling_factor,
                bench_options,
            )
        end
    end
end
end
