"""
    Definitions

Constants and parameter types shared by the whole package: tree configurations, flow
directions, ODE group codes and the parameter structs passed as `p` to the ODE functions
in pharmacokinetic_models.jl.
"""
module Definitions

    export TreeDefinitions,
        FlowDirections, 
        TreeStructure,  
        ODE_groups, 
        vascular_tree_parameters,
        terminal_parameters

    using Parameters

    using ..Paths: 
        PROJECT_ROOT, 
        JULIA_RESULTS_DIR

    """
    TreeDefinitions()

    Supported tree configurations. `vascular_trees` maps a configuration name to its trees,
    grouped as `:inflow_trees` and `:outflow_trees`:

    | Configuration      | Inflow trees | Outflow trees |
    |:-------------------|:-------------|:--------------|
    | `"Rectangle_trio"` | `A`, `P`     | `V`           |
    | `"Rectangle_quad"` | `A`, `P`     | `V`, `B`      |

    `trees.vascular_trees[tree_configuration]` is passed as `tree_components` to [`TreeStructure`](@ref).
    """
    @with_kw struct TreeDefinitions
        vascular_trees::Dict{String,Dict{Symbol,Vector{String}}} = Dict(
            "Rectangle_trio" => Dict(:inflow_trees => ["A", "P"], :outflow_trees => ["V"]),
            "Rectangle_quad" =>
                Dict(:inflow_trees => ["A", "P"], :outflow_trees => ["V", "B"]), #["P", "A", "V", "B"]
        )
    end

    """
    FlowDirections()

    Tree ids by direction of blood flow: `inflow_trees` (`A` arterial, `P` portal) carry
    blood into the liver and `outflow_trees` (`V` venous, `B` biliary) carry it out. Outflow
    trees have their edges reversed during processing (see `prepare_graph_info`).
    """
    @with_kw struct FlowDirections{NI,NO}
        inflow_trees::NTuple{NI,String} = ("A", "P")
        outflow_trees::NTuple{NO,String} = ("V", "B")
    end

    """
    TreeStructure(; tree_configuration, n_node)

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
    Base.@kwdef struct TreeStructure{S<:String}
        tree_configuration::S
        n_node::Int
        graph_id::S = "$(tree_configuration)_$(n_node)"
        tree_components::Dict{Symbol,Vector{S}}
        vascular_trees::Vector{S} = reduce(vcat, values(tree_components))
        GRAPH_DIR::S = joinpath(PROJECT_ROOT, JULIA_RESULTS_DIR, tree_configuration, graph_id)
    end

    """
    ODE_groups()

    Integer code of every edge, stored in the `ODE_groups` column. The ODE functions use it
    to choose the equation of the edge:

    | Field         | Code | Edge                                                         |
    |:--------------|:-----|:-------------------------------------------------------------|
    | `marginal`    | `0`  | `(0, root)`, connects the tree to the outside (dose input)   |
    | `other`       | `1`  | inner edges of the tree                                       |
    | `preterminal` | `2`  | edge that ends in a leaf node, connected to a terminal       |
    | `terminal`    | `3`  | self-edge `(n, n)` of a leaf node, the terminal compartment  |

    Stored as `Int16` to keep the `.arrow` tables and `p.ODE_groups` small.
    """
    @with_kw struct ODE_groups
        marginal::Int16 = 0
        preterminal::Int16 = 2
        terminal::Int16 = 3
        other::Int16 = 1
    end

    """
    vascular_tree_parameters(id, is_inflow, species_ids, flow_values, volume_values,
                             ODE_groups, pre_elements, post_elements)

    Parameters of one vessel tree (`p` for `jf_inflow!` / `jf_outflow!`), built by
    `get_tree_parameters` from `<tree>.arrow`. All vectors have one entry per edge, in the
    row order of the table.

    # Fields
    - `id`: tree id, e.g. `"A"`.
    - `is_inflow`: `true` for inflow trees, chooses `jf_inflow!` or `jf_outflow!`.
    - `species_ids`: names of the state variables, `<tree>_<source>_<target>`.
    - `flow_values`: flow of the edge [L/min], already multiplied by `flow_scaling_factor`.
    - `volume_values`: volume of the edge [L].
    - `ODE_groups`: group code of the edge, see [`ODE_groups`](@ref).
    - `pre_elements`, `post_elements`: indices of the predecessor / successor edges.
    """
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

    """
    terminal_parameters(; id, species_ids, flow_values, volumes)

    Parameters of the terminal nodes (`p` for `jf_terminal!`), built by `get_terminal_parameters`
    from `T.arrow`. Matrices have shape `(n_inflow+1) × n_terminals`: column k = terminal k,
    row 1 = the terminal itself, rows 2… = one row per inflow tree.

    # Fields
    - `id`: always `"T"`.
    - `species_ids`: names of the state variables; row 1 `T_k`, then the inflow species
    connected to terminal k.
    - `flow_values`: row 1 = summed outflow of the terminal, rows 2… = inflow from each inflow tree.
    - `volumes`: volume of one terminal (the same for all terminals).

    # Preallocated work arrays (derived, DO NOT PASS)
    - `terminal_matrix_size`: `size(species_ids)`.
    - `terminal_inflow`: `(n_inflow) × n_terminals`, inflow flow × concentration.
    - `terminal_outflow`: `1 × n_terminals`, outflow flow × concentration.
    - `terminal_difference`: `1 × n_terminals`, inflow - outflow of every terminal.
    """
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
end