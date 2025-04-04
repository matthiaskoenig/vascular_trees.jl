module Julia_from_jgraph
"""
Module which contains functions to load .arrow files (contain information of julia graph needed for ODE model)
    and to prepare vector of initial values and parameters for the ODE function.

Input: GRAPH_DIR (e to Rectangle_quad_10),
       vascular_tree (ex. "A") - id of a single tree
Output: u0 (vector of initial values) and 
        p (tuple that contains parameters for the differential equations)
"""

export get_ODE_parameters, get_initial_values

using DataFrames, Arrow
# include("../utils.jl")
using ..Utils.Definitions: tree_definitions, terminal_parameters, vascular_tree_parameters

trees = tree_definitions()

# using InteractiveUtils

function get_ODE_parameters(tree_info, graph_subsystem::String, flow_scaling_factor::AbstractFloat)
    # get path to the graph
    GRAPH_PATH::String = normpath(joinpath(tree_info.GRAPH_DIR, "graphs/$(graph_subsystem).arrow"))

    if graph_subsystem != "T"
        p = get_graph_parameters(GRAPH_PATH, flow_scaling_factor)
    else
        n_inflow::Integer = length(tree_info.tree_components[:inflow_trees])
        p = get_graph_parameters(GRAPH_PATH, n_inflow, flow_scaling_factor)
    end

    return p
end

function get_initial_values(species::Array{T}) where {T<:Number}
    u0 = zeros(size(species))
    return u0
end

function get_graph_parameters(GRAPH_PATH::String, flow_scaling_factor::AbstractFloat)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    p = vascular_tree_parameters(
        graph.vascular_tree_id[1],
        graph.is_inflow[1],
        Vector(graph.species_ids),
        Vector(graph.flows) .* flow_scaling_factor,
        Vector(graph.volumes),
        Vector(graph.ODE_groups),
        Vector(Vector.(graph.pre_elements)),
        Vector(Vector.(graph.post_elements)),
    )

    return p
end

function get_graph_parameters(GRAPH_PATH::String, n_inflow::Integer, flow_scaling_factor::AbstractFloat)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    p = terminal_parameters(; id = "T", species_ids = Array{String}(graph[1:(n_inflow+1), :]), flow_values = Array{Float64}(graph[(n_inflow+2):(n_inflow+2+n_inflow), :]) .* flow_scaling_factor, volumes = Float64((graph[end, 1])))

    return p
end

end
