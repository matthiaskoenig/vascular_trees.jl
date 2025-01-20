module Julia_from_jgraph
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Graph_structure: graph_frame
    import .Utils.Options: graph_options

    include("./julia_models.jl")
    import .Julia_models: jf_dxdt!, str_jf_dxdt

    using JLD2, OrdinaryDiffEq

    # ============ Specify options
    g_options::graph_options = graph_options(
        n_nodes=[10],  #750, 1000, 1250, 1500
        tree_ids=[
            # "Rectangle_single_inflow",
            # "Rectangle_quad",
            "Rectangle_trio",
            ],
    )

    function __init__()
        for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
            x0, p = get_ODE_components(tree_id=tree_id, n_node=n_node)
            x0_str::Vector{String} = p[7]
            dx_str, dx_vstr = str_jf_dxdt(["" for element_id in x0_str], ["" for element_id in x0_str], x0_str, p, 0.0)
            dx_str = ["d($(element_id))/dt = $(dx_str[ke])" for (ke, element_id) in enumerate(x0_str)]
            dx = jf_dxdt!([0.0 for _ in eachindex(x0)], x0, p, 0.0)

            printstyled("------------------------------------------------------------------------------------\n"; color = 124)
            printstyled("   Generating equations for arterial inflow from $tree_id, № of nodes = $(n_node)   \n"; color = 9)
            printstyled("------------------------------------------------------------------------------------\n"; color = 124)

            println("")
            printstyled("$(x0_str[length(x0_str)])\n"; color = 51)
            println("Initial value: $(x0[length(x0_str)])")
            printstyled("Equation: $(dx_str[length(x0_str)])"; color = :magenta)
            println("")
            printstyled("Equation: $(dx_vstr[length(x0_str)]) = $(dx[length(x0_str)])\n"; color = :blue)

            for i in eachindex(dx)
                if i != length(x0_str)
                    println("")
                    printstyled("$(x0_str[i])\n"; color = 51)
                    println("Initial value: $(x0[i])")
                    printstyled("Equation: $(dx_str[i]) = $(dx[i])\n"; color = :magenta)
                    printstyled("Equation: $(dx_vstr[i]) = $(dx[i])\n"; color = :blue)
                end
            end

        prob = ODEProblem(jf_dxdt!, 
            x0,
            (0.0, 10.0),
            p)

        sol = solve(
            prob, 
            Tsit5() # Rosenbrock23(), # Tsit5(), # CVODE_BDF
            )
        end
    end

    function get_ODE_components(; tree_id::String, n_node::Int32)::Tuple{Vector{Float64}, Tuple}
        graph = load_graph(tree_id, n_node)
        p::Tuple{Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Float64}, Vector{Float64}, Vector{String}, Vector{String}, Vector{String}} = (
            graph.edges,
            graph.terminal_edges,
            graph.start_edge,
            graph.preterminal_edges,
            graph.flows,
            graph.volumes,
            graph.element_ids,
            graph.flow_ids,
            graph.volume_ids
        )
        x0::Vector{Float64} = zeros(length(p[1]))
        set_initial_values!(x0, 1.0)
        return x0, p
    end

    function load_graph(tree_id::String, n_node::Int32)::Any
        # get graph id for correct path definition
        graph_id::String = "$(tree_id)_$(n_node)"
        # get path to the graph
        GRAPH_PATH::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs/A.jld2"))
        graph_file::JLD2.JLDFile{JLD2.MmapIO} = jldopen(GRAPH_PATH,"r"; typemap=Dict("Main.Process_julia_graph.Utils.Graph_structure.graph_frame" => graph_frame))
        graph = graph_file["graph"]

        return graph
    end

    function set_initial_values!(x0::Vector{Float64}, initial_value::Float64)
        x0[length(x0)] = initial_value
    end

end