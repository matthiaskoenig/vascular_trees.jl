module Julia_from_jgraph
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Definitions: tree_definitions, graph_frame
    import .Utils.Options: graph_options

    include("./julia_models.jl")
    import .Julia_models: jf_dxdt!

    using JLD2, OrdinaryDiffEq, Plots, DiffEqCallbacks
    using TimerOutputs
    using Sundials

    # ============ Specify options
    g_options::graph_options = graph_options(
        n_nodes=[1000],  #750, 1000, 1250, 1500
        tree_ids=[
            "Rectangle_quad",
            # "Rectangle_trio",
            ],
    )

    # Already specified in utils.jl
    trees::tree_definitions = tree_definitions()

    function __init__()
        for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
            to = TimerOutput()
            printstyled("------------------------------------------------------------------------------------\n"; color = 124)
            printstyled("   $(tree_id), № of nodes = $(n_node)   \n"; color = 9)
            printstyled("------------------------------------------------------------------------------------\n"; color = 124)
            vessel_tree = "A"
            x0, p = get_ODE_components(tree_id, n_node, vessel_tree)
            prob = ODEProblem(jf_dxdt!, 
                x0,
                (0.0, 10.0/60.0),
                p)
            @timeit to "1" sol = solve(
                prob, 
                Tsit5(),
                # CVODE_BDF(),
                # callback=cb_variant2
                abstol=1e-6,
                reltol=1e-6 # Rosenbrock23(), # Tsit5(), # CVODE_BDF
                )

            show(to, sortby=:firstexec)
            display(plot(sol))
        end
    end

    function get_ODE_components(tree_id::String, n_node::Int32, vessel_tree::String)::Tuple{Vector{Float64}, Tuple}
        # get graph id for correct path definition
        graph_id::String = "$(tree_id)_$(n_node)"
        # get path to the graph
        GRAPH_PATH::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs/$(vessel_tree).jld2"))

        printstyled("   Generating ODE components for $(vessel_tree)   \n"; color = 130)
        graph = load_graph(GRAPH_PATH)
        p = (
            graph.vascular_tree_id,
            graph.is_inflow,
            graph.edges,
            graph.terminal_edges,
            # graph.start_edge,
            #graph.preterminal_edges,
            graph.flows,
            graph.volumes,
            # graph.element_ids,
            # graph.flow_ids,
            # graph.volume_ids
            graph.group,
            graph.pre_elements,
            graph.post_elements
        )
        x0::Vector{Float64} = zeros(length(p[3]))
        (p[1] == "A") && (set_initial_values!(x0, 1.0))
        # set_initial_values!(x0, 1.0, p)
        return x0, p
    end

    function load_graph(GRAPH_PATH::String)::Any
        graph_file::JLD2.JLDFile{JLD2.MmapIO} = jldopen(GRAPH_PATH,"r"; typemap=Dict("Main.Process_julia_graph.Utils.Definitions.graph_frame" => graph_frame))
        graph = graph_file["graph"]

        return graph
    end

    function set_initial_values!(x0::Vector{Float64}, initial_value::Float64) # p::Tuple{Bool, Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Float64}, Vector{Float64}}
        # if p[1]
        #     x0[length(x0)] = initial_value
        # else
        #     x0[length(x0)-length(p[3]):length(x0)-1] = [initial_value for _ in eachindex(p[3])]
        # end
        x0[length(x0)] = initial_value
    end

end