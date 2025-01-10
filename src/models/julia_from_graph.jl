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
    1. x0, p

    FIXME matrices are two large. Once fully created matrix A should be transformed to sparse adjacency matrix.
    """

    # https://juliagraphs.org/Graphs.jl/dev/
    using CSV, DataFrames, Graphs, EzXML, ParameterizedFunctions, GraphDataFrameBridge, MetaGraphs, OrdinaryDiffEq, TimerOutputs, BenchmarkTools, SparseArrays
    # using GraphPlot
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Options: graph_options, edge_options
    include("./julia_models.jl")
    import .Julia_models: f_dxdt!

    export get_ODE_components

    function __init__()
        x0, p = get_ODE_components(tree_id="Rectangle_quad", n_node=Int32(100))
        dx::Matrix{Float64} = zeros(size(@view p[:, 1:Int(size(p, 2)/4)]))
        f_dxdt!(dx, x0, p, 0.0)
    end

    function get_ODE_components(; tree_id::String, n_node::Int32)::Tuple{Matrix{Float64}, Matrix{Float64}}
        # to = TimerOutput()
        edges_df, graph = read_graph(tree_id=tree_id, n_node=n_node)
        p = collect_graph_characteristics(edges_df=edges_df, graph=graph)
        x0::Vector{Float64} = zeros(size(p, 1))
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

        # calculation of terminal volume
        volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000  # [cm^3] -> [l]
        n_terminals::Int32 = nrow(edges_df[edges_df.terminal, :])
        volume_terminal::Float64 = volume_geometry / n_terminals

        A::Matrix{Float64} = convert(Matrix{Float64}, (adjacency_matrix(graph)))
        modify_adjacency_matrix!(A, graph)
        
        edges::Vector{CartesianIndex{2}} = findall(!iszero, A)
        p::Matrix{Float64} = zeros(length(edges), 5)
        collect_parameters_values!(p, edges, graph, volume_terminal)

        @show p

        return p
    end

    function modify_adjacency_matrix!(A::Matrix{Float64}, graph::MetaDiGraph{Int64, Float64})
        @inbounds for (index, value) in pairs(IndexCartesian(), A)
            source_id, target_id = Tuple.(index)
            ((target_id == source_id) && (startswith(props(graph, target_id)[:name], "T_"))) && (value = 1.0)
        end
    end

    function collect_parameters_values!(p::Matrix{Float64}, edges::Vector{CartesianIndex{2}}, graph::MetaDiGraph{Int64, Float64}, volume_terminal::Float64)
        # create a parameter vector, which will contain:
        # esges (source id, target id), volume values, flow_values, whether edge is from inflow system or not
        # THIS ORDER IS CRUCIAL
        sources = @view p[:, 1]
        targets = @view p[:, 2]
        volume_values = @view p[:, 3]
        flow_values = @view p[:, 4]
        is_inflow = @view p[:, 5]
        @inbounds for (ke, edge) ∈ enumerate(edges)
            source_id, target_id = Tuple.(edge)
            sources[ke] = source_id
            targets[ke] = target_id
            if target_id == source_id
                volume_values[ke] = volume_terminal
                is_inflow[ke] = 0.0
            else 
                volume_values[ke] = π * props(graph, source_id, target_id)[:radius]^2 * props(graph, source_id, target_id)[:length]
                flow_values[ke] = props(graph, source_id, target_id)[:flow]
                is_inflow[ke] = Float64(props(graph, source_id, target_id)[:is_inflow])
            end
        end
    end

end
