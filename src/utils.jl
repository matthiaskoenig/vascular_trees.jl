module Utils

    RESULTS_DIR = "results"
    JULIA_RESULTS_DIR = RESULTS_DIR * "/julia_vessel_trees"

    module Options
        export graph_options, simulations_options, benchmark_options

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
            plot_running_times::Bool = false
            n_iterations::Int16 = 5
        end
        
    end

    module Benchmarking
        export save_times_as_csv ,plot_running_times

        using CSV
        using TimerOutputs
        using DataFrames
        using PlotlyJS
        import ..JULIA_RESULTS_DIR

        function save_times_as_csv(; times::TimerOutput, n_species)
            individual_times = TimerOutputs.todict(times)["inner_timers"]
            labels::Array{String} = []
            n_calls::Array{Int16} = []
            times_min::Array{Float64} = [] 
            n_sp::Array{Int32} = []
            for graph_id ∈ keys(individual_times)
                n_call = parse(Int16, split(string(graph_id), "_")[4])
                time_min = get(get(individual_times, graph_id, NaN), "time_ns", NaN) / 60 * 10^-9
                push!(n_calls, n_call)
                push!(times_min, time_min)
                push!(labels, chopsuffix(string(graph_id), "_$(n_calls)"))
                push!(n_sp, n_species[n_call])
            end

            df = DataFrame([:graph_id=>labels, :n_calls=>n_calls, :times_min=>times_min, :n_species=>n_sp])

            CSV.write(joinpath(JULIA_RESULTS_DIR, "running_times.csv"), df)
        end

        function plot_running_times()

            running_times = CSV.read(joinpath(JULIA_RESULTS_DIR, "running_times.csv"))

            plt = plot(
                scatter(running_times, x=:n_calls, y=:n_species, z=:times_min, 
                        color=:graph_id, mode="markers",
                        type="scatter3d", labels=:graph_id),
                Layout(showlegend=true)
            )
            display(plt)

        end
    end
end