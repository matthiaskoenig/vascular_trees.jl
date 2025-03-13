module Julia_from_jgraph
"""
Module which contains functions to load .arrow files (contain information of julia graph needed for ODE model)
    and to prepare vector of initial values and parameters for the ODE function.

Input: tree_id (ex. "Rectangle_quad"), 
       n_node (ex. 10),
       vascular_tree (ex. "A") - id of a single tree
Output: u0 (vector of initial values) and 
        graph_p (structure that contains parameters for the differential equations)

This module is not used by itself. It's functions are imported and used by the functions in simulation_helpers.jl.
To use this module itsef for simulations - uncomment __init__ function (will be deleted in the future).
"""

export get_ODE_components

include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Definitions: tree_definitions
import .Utils.Options: graph_options

# include("./julia_models.jl")
# import .Julia_models: jf_dxdt!

using DataFrames, OrdinaryDiffEq, Plots, DiffEqCallbacks
using TimerOutputs, InteractiveUtils
using Sundials, Arrow

# ============ Specify options
g_options::graph_options = graph_options(
    n_nodes = [10],  #750, 1000, 1250, 1500
    tree_ids = [
        "Rectangle_quad",
        # "Rectangle_trio",
    ],
)

# Already specified in utils.jl
trees::tree_definitions = tree_definitions()

# function __init__()
#     for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
#         to = TimerOutput()
#         vascular_tree = "A"
#         #@code_warntype get_graph_components(tree_id, n_node, vascular_tree)
#         u0, graph_p = get_graph_components(tree_id, n_node, vascular_tree)
#         p = (graph_p.is_inflow, graph_p.flows, graph_p.volumes, graph_p.ODE_groups, graph_p.pre_elements, graph_p.post_elements)
#         #jf_dxdt!([0.0 for _ in eachindex(u0)], u0, p, 10.0)
#         prob = ODEProblem(jf_dxdt!, u0, (0.0, 10.0 / 60.0), p)
#         sol = solve(
#             prob,
#             Tsit5(),
#             # CVODE_BDF(),
#             # callback=cb_variant2
#             abstol = 1e-8,
#             reltol = 1e-8, # Rosenbrock23(), # Tsit5(), # CVODE_BDF
#         )

#         # show(to, sortby = :firstexec)
#         # display(plot(sol))
#     end
# end

function get_ODE_components(tree_id::String, n_node::Integer, vascular_tree::String)
    # get graph id for correct path definition
    graph_id::String = "$(tree_id)_$(n_node)"
    # get path to the graph
    GRAPH_PATH::String = normpath(
        joinpath(
            @__FILE__,
            "../../..",
            JULIA_RESULTS_DIR,
            tree_id,
            graph_id,
            "graphs/$(vascular_tree).arrow",
        ),
    )

    if vascular_tree != "T"
        parameters = get_graph_parameters(GRAPH_PATH)
        u0 = zeros(length(parameters[3]))
    else
        n_inflows::Integer = length(trees.vascular_trees[tree_id][:inflow_trees])
        parameters = get_graph_parameters(GRAPH_PATH, n_inflows)
        u0 = zeros(size(parameters[3]))
    end

    return u0, parameters
end

function get_graph_parameters(GRAPH_PATH::String)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    graph_parameters = (
        graph.vascular_tree_id[1],
        graph.is_inflow[1],
        Vector(graph.element_ids),
        Vector(graph.flows),
        Vector(graph.volumes),
        Vector(graph.ODE_groups),
        Vector(graph.pre_elements),
        Vector(graph.post_elements)
    )

    return graph_parameters
end

function get_graph_parameters(GRAPH_PATH::String, n_inflows::Integer)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    graph_parameters = (
        "T",
        Array{String}(graph[1:(n_inflows+1), :]), #x_affiliations
        Array{AbstractFloat}(graph[(n_inflows+2):(n_inflows+2+n_inflows), :]), #terminal_nodes.flow_values,
        # Vector(graph[(n_inflows+4):(n_inflows+5), 1]), #terminal_nodes.flow_affiliations,
        AbstractFloat((graph[end, 1]))  #terminal_nodes.volumes
    )

    return graph_parameters
end

end
