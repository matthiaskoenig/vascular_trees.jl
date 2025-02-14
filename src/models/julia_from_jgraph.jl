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

# function __init__()
#     for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
#         to = TimerOutput()
#         vessel_tree = "A"
#         # @code_warntype get_ODE_components(tree_id, n_node, vessel_tree)
#         x0, p = get_ODE_components(tree_id, n_node, vessel_tree)
#         # jf_dxdt!([0.0 for _ in eachindex(x0)], x0, p, (0.0, 10.0 / 60.0))
#         # prob = ODEProblem(jf_dxdt!, x0, (0.0, 10.0 / 60.0), p)
#         # @timeit to "1" sol = solve(
#         #     prob,
#         #     Tsit5(),
#         #     # CVODE_BDF(),
#         #     # callback=cb_variant2
#         #     abstol = 1e-8,
#         #     reltol = 1e-8, # Rosenbrock23(), # Tsit5(), # CVODE_BDF
#         # )

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

    graph = load_graph(GRAPH_PATH)
    p = get_parameters_value(graph)

    x0::Vector{Float64} = zeros(length(p[3]))
    (p[1] == "A") && (set_initial_values!(x0, 1.0))
    return x0, p
end

function load_graph(GRAPH_PATH::String)::DataFrame
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    select!(
        graph,
        [
            :vascular_tree_id,
            :is_inflow,
            :flows,
            :volumes,
            :ODE_groups,
            :pre_elements,
            :post_elements,
        ],
    )

    return graph
end

function set_initial_values!(x0::Vector{Float64}, initial_value::Float64)
    x0[length(x0)] = initial_value
end

function get_parameters_value(graph::DataFrame)
    vascular_tree_id::String = graph.vascular_tree_id[1]
    is_inflow::Bool = graph.is_inflow[1]
    flows::Vector{Float64} = graph.flows
    volumes::Vector{Float64} = graph.volumes
    ODE_groups::Vector{Int32} = graph.ODE_groups
    pre_elements::Vector{Vector{Int32}} = graph.pre_elements
    post_elements::Vector{Vector{Int32}} = graph.post_elements

    p = (
        vascular_tree_id,
        is_inflow,
        flows,
        volumes,
        ODE_groups,
        pre_elements,
        post_elements,
    )
    return p
end

end
