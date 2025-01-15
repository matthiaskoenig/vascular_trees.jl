module Julia_from_jgraph
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Graph_structure: graph_frame
    import .Utils.Options: graph_options

    include("./julia_models.jl")
    import .Julia_models: jf_dxdt!

    using JLD2

    # ============ Specify options
    g_options::graph_options = graph_options(
        n_nodes=[10],  #750, 1000, 1250, 1500
        tree_ids=[
            # "Rectangle_quad",
            "Rectangle_trio",
            ],
    )

    function __init__()
        for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
            x0, p = get_ODE_components(tree_id=tree_id, n_node=n_node)
            jf_dxdt!([0.0 for _ in 1:length(p[1])], x0, p, 0.0)
        end
    end

    function get_ODE_components(; tree_id::String, n_node::Int32)::Tuple{Vector{Float64}, Tuple}
        graph = load_graph(tree_id, n_node)
        p::Tuple{Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Float64}, Vector{Float64}} = (
            graph.edges,
            graph.terminal_edges,
            graph.start_edge,
            graph.preterminal_edges,
            graph.flows,
            graph.volumes
        )
        x0::Vector{Float64} = zeros(length(p[1]))
        set_initial_values!(x0, 10.0, p[1], p[3])
        return x0, p
    end

    function load_graph(tree_id::String, n_node::Int32)::Any
        # get graph id for correct path definition
        graph_id::String = "$(tree_id)_$(n_node)"
        # get path to the graph
        GRAPH_PATH::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs/A.jld2"))
        graph_file::JLD2.JLDFile{JLD2.MmapIO} = jldopen(GRAPH_PATH,"r"; typemap=Dict("Main.A" => graph_frame))
        graph = graph_file["graph"]

        return graph
    end

    function set_initial_values!(x0::Vector{Float64}, initial_value::Float64, edges::Vector{Tuple{Int32, Int32}}, start_edge::Vector{Tuple{Int32, Int32}})
        start_edge_idx = findfirst(x -> x==start_edge[1], edges)
        x0[start_edge_idx] = initial_value
    end

end