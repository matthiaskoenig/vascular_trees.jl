module Definitions
    export tree_definitions, Tree_structure, flow_directions, ODE_groups, terminal_parameters, vascular_tree_parameters

    using Parameters
    using Revise
    using ..Paths: PROJECT_ROOT, JULIA_RESULTS_DIR

    @with_kw struct tree_definitions
        vascular_trees::Dict{String,Dict{Symbol,Vector{String}}} = Dict(
            "Rectangle_trio" => Dict(:inflow_trees => ["A", "P"], :outflow_trees => ["V"]),
            "Rectangle_quad" =>
                Dict(:inflow_trees => ["A", "P"], :outflow_trees => ["V", "B"]), #["P", "A", "V", "B"]
        )
    end

    """
    Tree_structure(; tree_configuration, n_node)

    Basic information about the tree that differs between its types (Rectangle_quad, trio, etc.)
    and which is used repeatedly in simulations. Only `tree_configuration` and `n_node` are meant
    to be given, the other fields are derived from them (DO NOT CHANGE).

    # Fields
    - `tree_configuration`: type of the tree, e.g. `"Rectangle_quad"`.
    - `n_node`: number of nodes of the graph.
    - `graph_id`: `"<tree_configuration>_<n_node>"`, name of the graph folder.
    - `tree_components`: tree ids of the configuration, grouped as `:inflow_trees` / `:outflow_trees`.
    - `vascular_trees`: flat list of all tree ids, e.g. `["A", "P", "V", "B"]`.
    - `GRAPH_DIR`: directory with the graph files; simulation results go to its `simulations` subfolder.
    """
    Base.@kwdef struct Tree_structure{S<:String}
        tree_configuration::S
        n_node::Int
        graph_id::S = "$(tree_configuration)_$(n_node)"
        tree_components::Dict{Symbol,Vector{S}}
        vascular_trees::Vector{S} = reduce(vcat, values(tree_components))
        GRAPH_DIR::S = joinpath(PROJECT_ROOT, JULIA_RESULTS_DIR, tree_configuration, graph_id)
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

    Base.@kwdef struct terminal_parameters{T<:AbstractFloat, I<:Int}
        id::String
        species_ids::Array{String, 2}
        flow_values::Array{T, 2}
        volumes::T
        terminal_matrix_size::Tuple{I, I} = size(species_ids)
        terminal_inflow::Array{T, 2} = zeros(terminal_matrix_size[1]-1, terminal_matrix_size[2])
        terminal_outflow::Array{T, 2} = zeros(1, terminal_matrix_size[2])
        terminal_difference::Array{T, 2} = zeros(1, terminal_matrix_size[2])
    end

    struct vascular_tree_parameters{T<:AbstractFloat, I<:Int}
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