module Benchmarking
    export save_times_as_csv

    using CSV
    using TimerOutputs
    using Tables
    
    using ..Paths: JULIA_RESULTS_DIR, BENCHMARKING_RESULTS_PATH

    function save_times_as_csv(
        times::TimerOutput,
        n_node::Int64,
        tree_configuration::String,
        n_species::Int64,
        solver_name::String,
        n_term::Int64,
    )

        total_times = TimerOutputs.todict(times)["inner_timers"]
        n_calls::Vector{Int} = []
        call_orders::Vector{String} = []
        times_ns::Vector{Float64} = []
        allocated_bytes::Vector{Float64} = []
        for graph_id ∈ keys(total_times)
            indiv_times = get(get(total_times, graph_id, NaN), "inner_timers", NaN)
            len = length(keys(indiv_times))

            n_call = get(get(total_times, graph_id, NaN), "n_calls", NaN)
            time_ns = get(get(total_times, graph_id, NaN), "time_ns", NaN)
            allocated_mem = get(get(total_times, graph_id, NaN), "allocated_bytes", NaN)
            push!(n_calls, n_call)
            push!(call_orders, "Total time (all iterations)")
            push!(times_ns, time_ns)
            push!(allocated_bytes, allocated_mem)
            for subgraph_id ∈ keys(indiv_times)
                n_call_minor = get(get(indiv_times, subgraph_id, NaN), "n_calls", NaN)
                time_ns = get(get(indiv_times, subgraph_id, NaN), "time_ns", NaN)
                allocated_mem = get(get(indiv_times, subgraph_id, NaN), "allocated_bytes", NaN)
                push!(n_calls, n_call_minor)
                push!(call_orders, "Call")
                push!(times_ns, time_ns)
                push!(allocated_bytes, allocated_mem)
            end
        end
        tree_configurations::Vector{String} =
            ["$(tree_configuration)" for _ ∈ eachindex(allocated_bytes)]
        n_nodes::Vector{Int} = [n_node for _ ∈ eachindex(allocated_bytes)]
        n_sp::Vector{Int} = [n_species for _ ∈ eachindex(allocated_bytes)]
        solver_names::Vector{String} = [solver_name for _ ∈ eachindex(allocated_bytes)]
        n_terminals::Vector{Int} = [n_term for _ ∈ eachindex(allocated_bytes)]
        table = (
            n_calls = n_calls,
            call_orders = call_orders,
            times_min = times_ns / 60 * 10^-9,
            n_species = n_sp,
            allocated_gbytes = allocated_bytes / 2^30,
            solver_names = solver_names,
            n_node = n_nodes,
            tree_configuration = tree_configurations,
            n_terminals = n_terminals,
        )

        if isfile(BENCHMARKING_RESULTS_PATH)
            kwargs = (writeheader = false, append = true, sep = ',')
        else
            kwargs = (writeheader = true, sep = ',')
        end
        CSV.write(BENCHMARKING_RESULTS_PATH, table; kwargs...)
    end

end