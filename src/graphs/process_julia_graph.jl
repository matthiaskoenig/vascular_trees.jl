module Process_julia_graph
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Graph_structure: graph_frame
    import .Utils.Options: graph_options

    include("./processing_helpers.jl")
    import .Processing_helpers: read_edges, read_nodes_attributes

    using CSV, DataFrames, DataFramesMeta

    #graph::graph_frame = graph_frame()

    function __init__()
        graph_structure, nodes_attrib = read_graph(tree_id="Rectangle_trio", n_node=Int32(10))
        create_graph_structure(graph_structure, nodes_attrib)
    end

    function read_graph(; tree_id::String, n_node::Int32)::Tuple{DataFrame, DataFrame}
        # get graph id for correct path definition
        graph_id::String = "$(tree_id)_$(n_node)"
        # get directory with graph's files
        GRAPH_DIR::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "julia"))
        
        # read, prepare and join files with information about edges and their attributes
        graph_structure::DataFrame = read_edges(GRAPH_DIR)
        # read and prepare csv file with nodes attributes
        nodes_attrib::DataFrame = read_nodes_attributes(GRAPH_DIR)

        return graph_structure, nodes_attrib
    end

    function create_graph_structure(graph_structure::DataFrame, nodes_attrib::DataFrame)
        # nodes::Vector{Int32}
        # nodes_coordinates::Vector{Tuple{Float64, Float64, Float64}}
        # edges::Vector{Tuple{Int32, Int32}}
        # terminal_edges::Vector{Tuple{Int32, Int32}}
        # start_edges::Tuple{Int32, Int32}
        # preterminal_edges::Vector{Tuple{Int32, Int32}}
        # flows::Vector{Float64}
        # volumes::Vector{Float64}
        terminal_edges = subset(graph_skeleton, :terminal => x -> x .== true, view=true)
        start_edge = subset(graph_skeleton, :start => x -> x .== true, view=true)
        preterminal_edges = subset(graph_skeleton, :preterminal => x -> x .== true, view=true)
        @show graph::graph_frame = graph_frame(
            nodes_attrib.ids,
            Tuple.(Tables.namedtupleiterator(nodes_attrib[:, [:x, :y, :z]])),
            Tuple.(Tables.namedtupleiterator(graph_structure[:, [:source_id, :target_id]])),
            Tuple.(Tables.namedtupleiterator(terminal_edges[:, [:source_id, :target_id]])),
            # Tuple.(Tables.namedtupleiterator(start_edge[:, [:source_id, :target_id]])),
            Tuple.(Tables.namedtupleiterator(preterminal_edges[:, [:source_id, :target_id]])),
            graph_structure.flow,
            graph_structure.volumes
        )
    end

end