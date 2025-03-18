module Julia_from_jgraph
"""
Module which contains functions to load .arrow files (contain information of julia graph needed for ODE model)
    and to prepare vector of initial values and parameters for the ODE function.

Input: GRAPH_DIR (e to Rectangle_quad_10),
       vascular_tree (ex. "A") - id of a single tree
Output: u0 (vector of initial values) and 
        p (tuple that contains parameters for the differential equations)
"""

export get_ODE_components

using DataFrames, Arrow
include("../utils.jl")
import .Utils.Definitions: tree_definitions

trees = tree_definitions()

# using InteractiveUtils

function get_ODE_components(tree_info, graph_subsystem::String)
    # get path to the graph
    GRAPH_PATH::String = normpath(joinpath(tree_info.GRAPH_DIR, "graphs/$(graph_subsystem).arrow"))

    if graph_subsystem != "T"
        p = get_graph_parameters(GRAPH_PATH)
    else
        n_inflow::Integer = length(tree_info.tree_components[:inflow_trees])
        p = get_graph_parameters(GRAPH_PATH, n_inflow)
    end
    
    u0 = zeros(size(p[3]))

    return u0, p
end

function get_graph_parameters(GRAPH_PATH::String)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    p = (
        graph.vascular_tree_id[1],
        graph.is_inflow[1],
        Vector(graph.species_ids),
        Vector(graph.flows),
        Vector(graph.volumes),
        Vector(graph.ODE_groups),
        Vector(graph.pre_elements),
        Vector(graph.post_elements),
    )

    return p
end

function get_graph_parameters(GRAPH_PATH::String, n_inflow::Integer)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    p = (
        "T",
        Array{String}(graph[1:(n_inflow+1), :]), #x_affiliations
        Array{AbstractFloat}(graph[(n_inflow+2):(n_inflow+2+n_inflow), :]), #terminal_nodes.flow_values,
        # Vector(graph[(n_inflow+4):(n_inflow+5), 1]), #terminal_nodes.flow_affiliations,
        AbstractFloat((graph[end, 1])),  #terminal_nodes.volumes
    )

    return p
end

end
