module Utils
export JULIA_RESULTS_DIR, MODEL_PATH

RESULTS_DIR::String = "results"
JULIA_RESULTS_DIR::String = RESULTS_DIR * "/julia_vessel_trees"
BENCHMARKING_RESULTS_PATH::String = joinpath(JULIA_RESULTS_DIR, "jrunning_times.csv")
MODEL_PATH::String = "src/models/pharmacokinetic_models.jl"


module Definitions
export tree_definitions, flow_directions, ODE_groups, terminal_parameters, vascular_tree_parameters

using Parameters
using Revise

@with_kw struct tree_definitions
    vascular_trees::Dict{String,Dict{Symbol,Vector{String}}} = Dict(
        "Rectangle_trio" => Dict(:inflow_trees => ["A", "P"], :outflow_trees => ["V"]),
        "Rectangle_quad" =>
            Dict(:inflow_trees => ["A", "P"], :outflow_trees => ["V", "B"]), #["P", "A", "V", "B"]
    )
end

@with_kw struct flow_directions
    inflow_trees::Tuple{String,String} = ("A", "P")
    outflow_trees::Tuple{String,String} = ("V", "B")
end

@with_kw struct ODE_groups
    marginal::Int16 = 0
    preterminal::Int16 = 2
    terminal::Int16 = 3
    other::Int16 = 1
end

Base.@kwdef struct terminal_parameters{T<:AbstractFloat, I<:Integer}
    id::String
    species_ids::Array{String}
    flow_values::Array{T}
    volumes::T
    terminal_matrix_size::Tuple{I, I} = size(species_ids)
    terminal_inflow::Array{T} = zeros(terminal_matrix_size[1]-1, terminal_matrix_size[2])
    terminal_outflow::Array{T} = zeros(1, terminal_matrix_size[2])
    terminal_difference::Array{T} = zeros(1, terminal_matrix_size[2])
end

struct vascular_tree_parameters{T<:AbstractFloat, I<:Integer}
    id::String
    is_inflow::Bool
    species_ids::Vector{String}
    flow_values::Vector{T}
    volume_values::Vector{T}
    ODE_groups::Vector{Int16}
    pre_elements::Vector{Vector{I}}
    post_elements::Vector{Vector{I}}
end

end

module Options
export tree_options, simulations_options, benchmark_options, solver_options, edge_options
using Parameters

using OrdinaryDiffEq
using Revise # this package must not be in final version
using Sundials

@with_kw struct tree_options
    n_nodes::Vector{Int64}
    tree_configurations::Vector{String}
end

Base.@kwdef struct simulations_options
    tspan::Tuple{Float64,Float64}
    steps::Integer
    dt::Float64 = tspan[2] / steps
    save_simulations::Bool
    benchmark::Bool
end

@with_kw struct benchmark_options
    save_running_times::Bool = false
    n_iterations::Int16 = 1
end

# https://docs.sciml.ai/DiffEqDocs/stable/solvers/split_ode_solve/
@with_kw struct solver_options
    solver = Tsit5()
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-8
    solver_name = "Tsit5" # "Tsit5"
end

end

module Benchmarking
export save_times_as_csv

using CSV
using TimerOutputs
using Tables
import ..JULIA_RESULTS_DIR, ..BENCHMARKING_RESULTS_PATH

function save_times_as_csv(
    times::TimerOutput,
    n_node::Int64,
    tree_configuration::String,
    n_species::Int64,
    solver_name::String,
    n_term::Int64,
)

    total_times = TimerOutputs.todict(times)["inner_timers"]
    n_calls::Vector{Integer} = []
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
    n_nodes::Vector{Integer} = [n_node for _ ∈ eachindex(allocated_bytes)]
    n_sp::Vector{Integer} = [n_species for _ ∈ eachindex(allocated_bytes)]
    solver_names::Vector{String} = [solver_name for _ ∈ eachindex(allocated_bytes)]
    n_terminals::Vector{Integer} = [n_term for _ ∈ eachindex(allocated_bytes)]
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

end
