module Utils

    RESULTS_DIR = "results"
    JULIA_RESULTS_DIR = RESULTS_DIR * "/julia_vessel_trees"
    BENCHMARKING_RESULTS_PATH = joinpath(JULIA_RESULTS_DIR, "running_times.csv")

    module Options
        export graph_options, simulations_options, benchmark_options, solver_options, edge_options, model_types

        using Parameters
        using OrdinaryDiffEq
        using Revise # this package must not be in final version

        @with_kw struct graph_options
            n_nodes::Vector{Int32}
            tree_ids::Vector{String}
            model_types::Vector{String}
        end

        @with_kw struct simulations_options
            tspan::Tuple{Float64, Float64}
            tpoints::AbstractRange = []
            save_simulations::Bool
            benchmark::Bool
        end

        @with_kw struct benchmark_options
            save_running_times::Bool = false
            n_iterations::Int16 = 5
        end

        # https://docs.sciml.ai/DiffEqDocs/stable/solvers/split_ode_solve/
        @with_kw struct solver_options
            solver = Tsit5()
            absolute_tolerance::Float64 = 1e-6
            relative_tolerance::Float64 = 1e-6
            dt::Float64 = 0.001
            solver_name:: String = "Tsit5"
        end

        @with_kw struct edge_options
            id::String
            source::String
            target::String
            volume_id::String
            volume_value::Float32
            flow_id::String
            flow_value::Float32
            terminal::Bool
            start::Bool
            inflow::Bool
            initial::Float32
        end

        @with_kw struct model_types
            templates::Vector{String} = [
                                        "vectorized!_cleaned_typed",
                                        "vectorized_typed",
                                        "vectorized", 
                                        "vectorized!", 
                                        "vectorized_cleaned", 
                                        "vectorized!_cleaned",
                                        "symbolic_MT"
                                        ]
            julia_model::Vector{String} = ["vectorized!_loop_typed"]
        end
        
    end

    module Benchmarking
        export save_times_as_csv

        using CSV
        using TimerOutputs
        using Tables
        import ..JULIA_RESULTS_DIR, ..BENCHMARKING_RESULTS_PATH

        function save_times_as_csv(; times::TimerOutput, n_species, model_type, solver_name)

            individual_times = TimerOutputs.todict(times)["inner_timers"]
            labels::Vector{String} = []
            graph_ids::Vector{String} = []
            n_calls::Vector{Int16} = []
            call_orders::Vector{String} = []
            times_ns::Vector{Float64} = [] 
            allocated_bytes:: Vector{Float64} = []
            n_sp::Vector{Int32} = []
            for graph_id ∈ keys(individual_times)
                without_compilation = get(get(individual_times, graph_id, NaN), "inner_timers", NaN)
                len = length(keys(without_compilation))

                n_call = get(get(individual_times, graph_id, NaN), "n_calls", NaN)
                time_ns = get(get(individual_times, graph_id, NaN), "time_ns", NaN) 
                allocated_mem = get(get(individual_times, graph_id, NaN), "allocated_bytes", NaN)
                push!(n_calls, n_call)
                push!(call_orders, "Time with f_dxdt precompilation, first call")
                push!(times_ns, time_ns)
                push!(allocated_bytes, allocated_mem)
                push!(labels, string(graph_id))
                push!(graph_ids, graph_id)
                push!(n_sp, n_species[1])
                for (ksb, subgraph_id) ∈ enumerate(keys(without_compilation))
                    n_call_minor = get(get(without_compilation, subgraph_id, NaN), "n_calls", NaN)
                    time_ns = get(get(without_compilation, subgraph_id, NaN), "time_ns", NaN)
                    allocated_mem = get(get(without_compilation, subgraph_id, NaN), "allocated_bytes", NaN)
                    push!(n_calls, n_call_minor)
                    if ksb == 1
                        push!(call_orders, "Time without f_dxdt precompilation, first call")
                    else
                        push!(call_orders, "Time without f_dxdt precompilation, call number > 1 (n=$(len-1))")
                    end
                    push!(times_ns, time_ns)
                    push!(allocated_bytes, allocated_mem)
                    push!(labels, string(subgraph_id))
                    push!(graph_ids, graph_id)
                    push!(n_sp, n_species[n_call])
                end
            end
            model_types = [model_type for _ ∈ 1:length(graph_ids)]
            solver_names = [solver_name for _ ∈ 1:length(graph_ids)]
            table = (graph_ids=graph_ids, 
                    labels=labels, 
                    n_calls=n_calls, 
                    call_orders=call_orders, 
                    times_min=times_ns/60*10^-9, 
                    n_species=n_sp, 
                    allocated_gbytes=allocated_bytes/10^9, 
                    model_types=model_types, 
                    solver_names=solver_names)

            if isfile(BENCHMARKING_RESULTS_PATH)
                kwargs = (writeheader = false, append = true, sep=',')
            else
                kwargs = (writeheader = true, sep=',')
            end
            CSV.write(BENCHMARKING_RESULTS_PATH, table; kwargs...)
        end

    end
end