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