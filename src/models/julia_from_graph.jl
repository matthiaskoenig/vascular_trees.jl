module Julia_from_graph
    """
    Module which contains functions for creating x0 and p vectors for ODESolve function.
        These matrices are needed for ODE functions from julia_models.jl.

    Also this module can be used for checking ODE functions from julia_models.jl.
        1. For this uncomment:
            #include("./julia_models.jl")
            #import .Julia_models: f_dxdt!
            __init__ function
        2. Specialize what graph do you want to use in  __init__ function


    Input:
    1. graph.csv - DataFrame with information about edges of the graph
       and their attributes

    Output:
    1. f_dxdt!

    FIXME matrices are two large. Once fully created matrix A should be transformed to sparse adjacency matrix.
    """

    # https://juliagraphs.org/Graphs.jl/dev/
    using CSV, DataFrames, Graphs, EzXML, ParameterizedFunctions, GraphDataFrameBridge, MetaGraphs, OrdinaryDiffEq, TimerOutputs, BenchmarkTools
    # using GraphPlot
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Options: graph_options, edge_options
    #include("./julia_models.jl")
    #import .Julia_models: f_dxdt!

    export get_ODE_components

    # function __init__()
    #     x0, p = get_ODE_components(tree_id="Rectangle_quad", n_node=Int32(100))
    #     dx::Matrix{Float64} = zeros(size(@view p[:, 1:Int(size(p, 2)/4)]))
    #     f_dxdt!(dx, x0, p, 0.0)
    # end

    function get_ODE_components(; tree_id::String, n_node::Int32)::Tuple{Matrix{Float64}, Matrix{Float64}}
        # to = TimerOutput()
        edges_df, graph = read_graph(tree_id=tree_id, n_node=n_node)
            p = collect_graph_characteristics(edges_df=edges_df, graph=graph)
        x0::Matrix{Float64} = zeros(size(@view p[:, 1:Int(size(p, 2)/4)]))
        return x0, p
    end

    function read_graph(; tree_id::String, n_node::Int32)::Tuple{DataFrame, MetaDiGraph{Int64, Float64}}
        # get graph id
        graph_id::String = "$(tree_id)_$(n_node)"
        # get path to the graph
        GRAPH_PATH::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs/graph.csv"))
        # read csv file with information about edges
        edges_df::DataFrame = DataFrame(CSV.File(GRAPH_PATH))
        # create graph from DataFrame
        graph::MetaDiGraph{Int64, Float64} = MetaDiGraph(edges_df, :source, :target,
                                                        edge_attributes=[:terminal, :start, :group, :is_inflow, :radius, :flow, :length])

        return edges_df, graph
    end

    function collect_graph_characteristics(; edges_df::DataFrame, graph::MetaDiGraph{Int64, Float64})::Matrix{Float64}

        # create a parameter vector, which will contain:
        # adjacency matrix, volume values, flow_values, whether edge is from inflow system or not
        # THIS ORDER IS CRUCIAL
        p::Array{Float64} = repeat(convert(Matrix{Float64}, (adjacency_matrix(graph))), outer=(1, 4))

        # create views for convenience
        # adjacency matrix for a graph, indexed by [u, v] vertices, so rows - sources, columns - targets
        A = @view p[:, 1:Int(size(p, 2)/4)]
        # volumes
        volume_values = @view p[:, size(A, 2)+1:size(A, 2)*2]
        # flows
        flow_values = @view p[:, size(A, 2)*2+1:size(A, 2)*3]
        # inflow system or not
        is_inflow = @view p[:, size(A, 2)*3+1:size(p, 2)]

        # calculation of terminal volume
        volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000  # [cm^3] -> [l]
        n_terminals::Int32 = nrow(edges_df[edges_df.terminal, :])
        volume_terminal::Float64 = volume_geometry / n_terminals


        @inbounds for (index, value) in pairs(IndexCartesian(), A)
            source_id, target_id = Tuple.(index)
            if value == 1.0
                volume_values[index] = π * props(graph, source_id, target_id)[:radius]^2 * props(graph, source_id, target_id)[:length]
                flow_values[index] = props(graph, source_id, target_id)[:flow]
                is_inflow[index] = Float64(props(graph, source_id, target_id)[:is_inflow])
            end
            if (target_id == source_id) && (startswith(props(graph, target_id)[:name], "T_"))
                value = 1.0
                volume_values[index] = volume_terminal
                is_inflow[index] = 0.0
            end
        end

        #@show p

        return p
    end

end
