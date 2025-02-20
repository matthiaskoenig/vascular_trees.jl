module Julia_from_jgraph
include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Definitions: tree_definitions
import .Utils.Options: graph_options

include("./julia_models.jl")
import .Julia_models: jf_dxdt!

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

struct to_collect
    vascular_tree_id::String
    is_inflow::Bool
    all_edges::Vector{Tuple{Int32, Int32}}
    flows::Vector{Float64}
    volumes::Vector{Float64}
    ODE_groups::Vector{Int32}
    pre_elements::Vector{Vector{Int32}}
    post_elements::Vector{Vector{Int32}}
end

# function __init__()
#     for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
#         to = TimerOutput()
#         vessel_tree = "A"
#         # @code_warntype get_ODE_components(tree_id, n_node, vessel_tree)
#         x0, graph_p = get_ODE_components(tree_id, n_node, vessel_tree)
#         p = (graph_p.is_inflow, graph_p.flows, graph_p.volumes, graph_p.ODE_groups, graph_p.pre_elements, graph_p.post_elements)
#         #jf_dxdt!([0.0 for _ in eachindex(x0)], x0, p, 10.0)
#         prob = ODEProblem(jf_dxdt!, x0, (0.0, 10.0 / 60.0), p)
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

function get_ODE_components(tree_id::String, n_node::Int32, vessel_tree::String)
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
            "graphs/$(vessel_tree).arrow",
        ),
    )

    graph_p = get_graph_parameters(GRAPH_PATH)

    x0::Vector{Float64} = zeros(length(graph_p.all_edges))
    # (graph_p.vascular_tree_id == "A") && (set_initial_values!(x0, 1.0))
    return x0, graph_p
end

function get_graph_parameters(GRAPH_PATH::String)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    graph_p = to_collect(
        graph.vascular_tree_id[1],
        graph.is_inflow[1],
        graph.all_edges,
        graph.flows,
        graph.volumes,
        graph.ODE_groups,
        graph.pre_elements,
        graph.post_elements,
    )

    return graph_p
end

function set_initial_values!(x0::Vector{Float64}, initial_value::Float64)
    x0[length(x0)] = initial_value
end


end
