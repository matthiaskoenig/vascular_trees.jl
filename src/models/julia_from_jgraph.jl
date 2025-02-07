module Julia_from_jgraph
include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Definitions: tree_definitions
import .Utils.Options: graph_options
import .Utils.Nice_printing: print_graph_info

include("./julia_models.jl")
import .Julia_models: jf_dxdt!

using DataFrames, OrdinaryDiffEq, Plots, DiffEqCallbacks
using TimerOutputs
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
#         x0, p = get_ODE_components(tree_id, n_node, vessel_tree)
#         prob = ODEProblem(jf_dxdt!, x0, (0.0, 10.0 / 60.0), p)
#         @timeit to "1" sol = solve(
#             prob,
#             Tsit5(),
#             # CVODE_BDF(),
#             # callback=cb_variant2
#             abstol = 1e-8,
#             reltol = 1e-8, # Rosenbrock23(), # Tsit5(), # CVODE_BDF
#         )

#         show(to, sortby = :firstexec)
#         #display(plot(sol))
#     end
# end

function get_ODE_components(
    tree_id::String,
    n_node::Int32,
    vessel_tree::String,
)::Tuple{Vector{Float64},Tuple}
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

    printstyled("   Generating ODE components for $(vessel_tree)   \n"; color = 130)
    graph = load_graph(GRAPH_PATH)
    p = (
        create_vector_from_column(graph.vascular_tree_id),
        create_vector_from_column(graph.is_inflow),
        create_vector_from_column(graph.all_edges),
        # graph.terminal_edges,
        # graph.start_edge,
        # graph.preterminal_edges,
        create_vector_from_column(graph.flows),
        create_vector_from_column(graph.volumes),
        # graph.element_ids,
        # graph.flow_ids,
        # graph.volume_ids
        create_vector_from_column(graph.ODE_groups),
        create_vector_from_column(graph.pre_elements),
        create_vector_from_column(graph.post_elements),
    )
    x0::Vector{Float64} = zeros(length(p[3]))
    (p[1][1] == "A") && (set_initial_values!(x0, 1.0))
    # set_initial_values!(x0, 1.0, p)
    return x0, p
end

function load_graph(GRAPH_PATH::String)::Arrow.Table
    graph = Arrow.Table(GRAPH_PATH)

    return graph
end

function set_initial_values!(x0::Vector{Float64}, initial_value::Float64)
    x0[length(x0)] = initial_value
end

function create_vector_from_column(column::Arrow.ArrowVector)::Vector
    return collect(skipmissing(column))
end

end
