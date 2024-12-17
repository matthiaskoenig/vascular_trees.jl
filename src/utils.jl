module Utils

    RESULTS_DIR = "results"
    JULIA_RESULTS_DIR = RESULTS_DIR * "/julia_vessel_trees"
    BENCHMARKING_RESULTS_PATH = joinpath(JULIA_RESULTS_DIR, "running_times.csv")

    module Options
        export graph_options, simulations_options, benchmark_options, solver_options

        using Parameters

        @with_kw struct graph_options
            n_nodes::Vector{Int32}
            tree_ids::Vector{String}
            file_suffix::String
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

        @with_kw mutable struct solver_options
            solver = Tsit5()
            absolute_tolerance::Float64 = 1e-6
            relative_tolerance::Float64 = 1e-6
        end
        
    end

    module Benchmarking
        export save_times_as_csv

        using CSV
        using TimerOutputs
        using Tables
        import ..JULIA_RESULTS_DIR, ..BENCHMARKING_RESULTS_PATH

        function save_times_as_csv(; times::TimerOutput, n_species)

            individual_times = TimerOutputs.todict(times)["inner_timers"]
            labels::Array{String} = []
            graph_ids::Array{String} = []
            n_calls::Array{Int16} = []
            call_orders::Array{String} = []
            tree_ids::Array{Int16} = []
            times_ns::Array{Float64} = [] 
            allocated_bytes:: Array{Float64} = []
            n_sp::Array{Int32} = []
            for graph_id ∈ keys(individual_times)
                without_compilation = get(get(individual_times, graph_id, NaN), "inner_timers", NaN)

                n_call = get(get(individual_times, graph_id, NaN), "n_calls", NaN)
                time_ns = get(get(individual_times, graph_id, NaN), "time_ns", NaN) 
                allocated_mem = get(get(individual_times, graph_id, NaN), "allocated_bytes", NaN)
                push!(n_calls, n_call)
                push!(call_orders, "With module inclusion (n=$(n_call / length(keys(without_compilation)))")
                push!(times_ns, time_ns)
                push!(allocated_bytes, allocated_mem)
                push!(labels, string(graph_id))
                push!(graph_ids, graph_id)
                push!(n_sp, n_species[1])
                push!(tree_ids, chopsuffix(graph_id, r"_\d"))
                for subgraph_id ∈ keys(without_compilation)
                    n_call_minor = get(get(without_compilation, subgraph_id, NaN), "n_calls", NaN)
                    time_ns = get(get(without_compilation, subgraph_id, NaN), "time_ns", NaN)
                    allocated_mem = get(get(without_compilation, subgraph_id, NaN), "allocated_bytes", NaN)
                    push!(n_calls, n_call_minor)
                    push!(call_orders, "Without module inclusion (n=$n_call)")
                    push!(times_ns, time_ns)
                    push!(allocated_bytes, allocated_mem)
                    push!(labels, string(subgraph_id))
                    push!(graph_ids, graph_id)
                    push!(n_sp, n_species[n_call])
                    push!(tree_ids, chopsuffix(graph_id, r"_\d"))
                end
            end
            table = (graph_ids=graph_ids, labels=labels, tree_ids=tree_ids, n_calls=n_calls, call_orders=call_orders, times_min=times_ns/60*10^-9, n_species=n_sp, allocated_gbytes=allocated_bytes/10^9)

            if isfile(BENCHMARKING_RESULTS_PATH)
                kwargs = (writeheader = false, append = true, sep=',')
            else
                kwargs = (writeheader = true, sep=',')
            end
            CSV.write(BENCHMARKING_RESULTS_PATH, table; kwargs...)
        end

    end
end