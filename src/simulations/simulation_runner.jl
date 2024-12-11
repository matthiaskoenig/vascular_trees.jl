
module Simulation_runner
    """
    Module that summarizes functions that load graph files created in Python (system of ODEs written in julia format, 
    initial values, parameters' values) and run simulations.

    Make sure that structure of vascular_trees.jl is as needed (somewhere it should be written)

    What should be specified in code below:
    1. g_options (graph options)
    2. sim_options (simulation options)

    Inputs:
    1. Graph files (ODE models written in julia format). Each graph may have model written as 
    in vectorized as in symbolic form.

    Outputs:
    1. simulations.csv - if sim_options.save_simulations = true
    2. running_times.csv - if sim_options.save_running_times = true

    TODO: run macros @timeit only when you need it - now it is so, but it is dirty
    FIXME: symbolic does not work - check template
    """

    include("../utils.jl")
    import .Utils.Options: graph_options, simulations_options

    include("simulation_helpers.jl")
    import .Simulation_helpers: create_simulations, create_benchmarked_simulations

    # ============ Specify options
    g_options::graph_options = graph_options(
        n_nodes=[10],
        tree_ids=[
            #"Rectangle_quad",
            "Rectangle_trio",],
        file_suffix="vectorized" #"vectorized" "symbolic"
    )

    tspan=(0.0, 10.0/60)
    sim_options::simulations_options = simulations_options(
        tspan=tspan,
        tpoints=range(tspan[1], stop=tspan[2], length=1001),
        save_simulations=true,
        benchmark=false
    )

    if sim_options.benchmark

        import .Utils.Options: benchmark_options

        bench_options:: benchmark_options = benchmark_options()
        create_benchmarked_simulations(g_options=g_options, sim_options=sim_options, bench_options=bench_options)
        # if bench_options.plot_running_times
        #     plot_running_times()
        # end
    else
        create_simulations(g_options=g_options, sim_options=sim_options)
    end

end

