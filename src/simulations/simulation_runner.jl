
module Simulation_runner
    """
    Module that summarizes functions that load graph files created in Python (system of ODEs written in julia format, 
    initial values, parameters' values) and runs simulations.

    Make sure that structure of vascular_trees.jl directory is as needed (somewhere it should be written)

    What must/may be specified in code below:
    1. must
        1.1. g_options (graph options)
        1.2. sim_options (simulation options)
    2. may (if you whant to geet additional data and do not like default variants)
        2.1. benchmark options 

    Inputs:
    1. Graph files (ODE models written in julia format). Each graph may have model written as 
    in different forms (see somewhere).

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
    import .Simulation_helpers: create_simulations, create_benchmarked_simulations

    # ============ Specify options
    g_options::graph_options = graph_options(
        n_nodes=[10, 30, 50, 100, 250, 500, 750, 1000, 1250, 1500, 1750],  # 10, 30, 50, 100, 250, 500, 750, 1000, 1250, 1500, 1750, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 200000, 300000, 400000, 500000, 1000000
        tree_ids=[
            "Rectangle_quad",
            # "Rectangle_trio",
            ],
        model_types=[
            # "loop!_python", # not working
            "loop!_julia",
            # "vectorized!_cleaned_typed",
            # "vectorized_typed",
            # "vectorized", 
            # "vectorized!", 
            # "vectorized_cleaned", 
            # "vectorized!_cleaned",
            # "symbolic_MT",
            ] 
    )

    tspan=(0.0, 10.0/59)
    sim_options::simulations_options = simulations_options(
        tspan=tspan,
        tpoints=range(tspan[1], stop=tspan[2], length=500),
        save_simulations=false,
        benchmark=true
    )

    sol_options:: solver_options = solver_options()

    if sim_options.benchmark

        import .Utils.Options: benchmark_options

        # do not write anything here in brackets if you are okay with default variant
        bench_options:: benchmark_options = benchmark_options(
            save_running_times=true
        )
        create_benchmarked_simulations(g_options=g_options, sim_options=sim_options, bench_options=bench_options, sol_options=sol_options)

    else
        create_simulations(g_options=g_options, sim_options=sim_options, sol_options=sol_options)
    end

end


