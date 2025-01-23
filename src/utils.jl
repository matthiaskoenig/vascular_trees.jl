module Utils

    RESULTS_DIR = "results"
    JULIA_RESULTS_DIR = RESULTS_DIR * "/julia_vessel_trees"
    BENCHMARKING_RESULTS_PATH = joinpath(JULIA_RESULTS_DIR, "jrunning_times.csv")

    
    module Definitions
        export tree_definitions, graph_frame

        using Parameters
        using Revise

        @with_kw struct tree_definitions
            vascular_trees::Dict{String, Vector{String}} = Dict(
                "Rectangle_trio" => ["A", "P", "V"],
                "Rectangle_quad" => ["A", "P", "V", "B"]
            )
            inflow_trees::Tuple{String, String} = ("A", "P")
            outflow_trees::Tuple{String, String} = ("V", "B")
        end
        
        struct graph_frame
            vascular_tree_id::String
            is_inflow::Bool
            nodes::Vector{Int32} # ids of all nodes
            nodes_coordinates::Vector{Tuple{Float64, Float64, Float64}} # position of each node (x, y, z)
            edges::Vector{Tuple{Int32, Int32}} # all edges
            terminal_edges::Vector{Tuple{Int32, Int32}} # edges between terminal nodes (the lowest nodes in network) and themselves
            start_edge::Vector{Tuple{Int32, Int32}} # MUST be changed to just Tuple # the highest node in network
            preterminal_edges::Vector{Tuple{Int32, Int32}} # edges between the terminal nodes and nodes on the level higher
            flows::Vector{Float64} # flow values
            volumes::Vector{Float64} # volume values
            element_ids::Vector{String}
            flow_ids::Vector{String}
            volume_ids::Vector{String}
        end

    end

    module Options
        export graph_options, simulations_options, benchmark_options, solver_options, edge_options, model_types
        using Parameters

        using OrdinaryDiffEq
        using Revise # this package must not be in final version

        @with_kw struct graph_options
            n_nodes::Vector{Int32}
            tree_ids::Vector{String}
            model_types::Union{Vector{String},Nothing} = nothing
        end

        @with_kw struct simulations_options
            tspan::Tuple{Float64, Float64}
            tpoints::AbstractRange = []
            save_simulations::Bool
            benchmark::Bool
        end

        @with_kw struct benchmark_options
            save_running_times::Bool = false
            n_iterations::Int16 = 2
        end

        # https://docs.sciml.ai/DiffEqDocs/stable/solvers/split_ode_solve/
        @with_kw struct solver_options
            solver = Tsit5()
            absolute_tolerance::Float64 = 1e-6
            relative_tolerance::Float64 = 1e-6
            dt::Float64 = 0.1
            solver_name:: String = "Tsit5" # "Tsit5"
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
            julia_model::Vector{String} = ["loop!_python", #not working
                                            "loop!_julia"]
        end
        
    end

    module Benchmarking
        export save_times_as_csv

        using CSV
        using TimerOutputs
        using Tables
        import ..JULIA_RESULTS_DIR, ..BENCHMARKING_RESULTS_PATH

        function save_times_as_csv(; times::TimerOutput, n_node::Int32, tree_id::String, n_species::Int64, model_type::String, solver_name::String, n_term::Int32)

            individual_times = TimerOutputs.todict(times)["inner_timers"]
            tree_ids::Vector{String} = []
            n_nodes::Vector{Int32} = []
            n_calls::Vector{Int16} = []
            call_orders::Vector{String} = []
            times_ns::Vector{Float64} = [] 
            allocated_bytes:: Vector{Float64} = []
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
                push!(tree_ids, tree_id)
                push!(n_nodes, n_node)
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
                    push!(tree_ids, tree_id)
                    push!(n_nodes, n_node)
                end
            end
            n_sp::Vector{Int64} = [n_species for _ ∈ 1:length(tree_ids)]
            model_types = [ "One_tree" for _ ∈ 1:length(tree_ids)] #
            solver_names = [solver_name for _ ∈ 1:length(tree_ids)]
            n_terminals::Vector{Int32} = [n_term for _ ∈ 1:length(tree_ids)]
            table = (n_calls=n_calls, 
                    call_orders=call_orders, 
                    times_min=times_ns/60*10^-9, 
                    n_species=n_sp, 
                    allocated_gbytes=allocated_bytes/10^9, 
                    model_types=model_types, 
                    solver_names=solver_names,
                    n_node=n_nodes,
                    tree_id=tree_ids,
                    n_terminals=n_terminals)

            if isfile(BENCHMARKING_RESULTS_PATH)
                kwargs = (writeheader = false, append = true, sep=',')
            else
                kwargs = (writeheader = true, sep=',')
            end
            CSV.write(BENCHMARKING_RESULTS_PATH, table; kwargs...)
        end

    end
end