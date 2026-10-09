"""
Load graph files in arrow format and run the coupled-tree simulations.

The simulations start as soon as this file is included
(e.g. `julia --project src/simulations/simulation_runner.jl`), using the options
defined in the options sections below.

Make sure that the structure of the vascular_trees.jl directory is as needed
(see "Inputs" for the expected graph location).

What must/may be specified in the code below:
1. must
    1.1. t_options (graph options)
    1.2. sim_options (simulation options)
2. may (if you want to add additional data and do not like default variants)
    2.1. benchmark options
    2.2. solver options
    2.3. additional solver options (these arguments are not mandatory)

Inputs:
1. Graph file (prepared information in arrow format (process_julia_graph.jl) from a graph
   generated in SyntheticVascularTrees.jl), expected in
   `results/julia_vessel_trees/<tree_configuration>/<tree_configuration>_<n_nodes>/`

Showed in the terminal:
1. Julia stuff (progress bar, @info messages)
2. benchmark results if sim_options.benchmark=true

Outputs:
1. `<subsystem>_simulations_<dt>_dt.csv` in the `simulations` folder of the graph directory,
   one file per vascular tree and one for the terminal part - if you do not do benchmarking
   and sim_options.save_simulations = true
2. jrunning_times.csv - if sim_options.benchmark and bench_options.save_running_times = true

TODO: run macros @timeit only when you need it - now it is so, but it is dirty
"""
module Simulation_Runner
    # === Imports ===
    using Revise
    using VascularTrees
    using VascularTrees.Helpers: run_simulations

    # Definitions of the available tree configurations (inflow / outflow trees)
    const trees::tree_definitions = tree_definitions()


    # ============================ Options ============================

    # === Graph options ===
    # options for graph, i.e., number of nodes and type of tree
    t_options = tree_options(
        n_nodes = [10, 100, 1000],  #1000, 10000, 100000, 100000
        tree_configurations = [
            "Rectangle_quad",
            "Rectangle_trio",
        ],
    )

    # factors by which the flow is scaled, 1.0 is the original flow
    flow_scaling_factors = [1.0] # 1.0/16, 1.0/8, 1.0/4, 1.0/2, 1.0, 1.0*2, 1.0*4, 1.0*8, 1.0*16

    # === Simulation options ===
    sim_options = simulations_options(
        tspan = (0.0, 16.0),  # [min]
        steps = 800,          # number of time steps, dt = tspan[2] / steps
        save_simulations = true,
        benchmark = true,
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
    bench_options = benchmark_options(save_running_times = false)
    

    # ============================== Run ==============================

    """
        simulate_all_trees()

    Run the simulation for every combination of tree configuration, number of nodes and flow
    scaling factor from the options above. Order: configurations (outermost), number of nodes,
    flow scaling factors (innermost).
    """
    function simulate_all_trees()
        # product varies its first iterator fastest, so the order is reversed here
        for (flow_scaling_factor, n_node, tree_configuration) ∈ Iterators.product(
            flow_scaling_factors,
            t_options.n_nodes,
            t_options.tree_configurations,
        )
            tree_info =
                Tree_structure(; tree_configuration = tree_configuration, n_node = n_node, tree_components = trees.vascular_trees[tree_configuration])
            @info "Working on $(tree_info.graph_id)"
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

    simulate_all_trees()

end
