module Process_julia_graph
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Graph_structure: graph_frame
    import .Utils.Options: graph_options

    include("./processing_helpers.jl")
    import .Processing_helpers: read_edges, read_nodes_attributes

    using DataFrames, JLD2, InteractiveUtils

    """
    TODO: Get rid of zeros in DataFrame
    """

    # ============ Specify options
    g_options::graph_options = graph_options(
        n_nodes=[10],  #750, 1000, 1250, 1500
        tree_ids=[
            "Rectangle_single_inflow",
            # "Rectangle_quad",
            # "Rectangle_trio",
            ],
    )

    function __init__()
        for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
            process_julia_graph(tree_id, n_node)
        end
    end

    function process_julia_graph(tree_id::String, n_node::Int32)
        # get graph id for correct path definition
        graph_id::String = "$(tree_id)_$(n_node)"
        graph_structure, nodes_attrib = read_graph(tree_id, graph_id)
        graph = create_graph_structure(graph_structure, nodes_attrib)
        save_graph!(graph, tree_id, graph_id)
    end

    function read_graph(tree_id::String, graph_id::String)::Tuple{DataFrame, DataFrame}
        # get directory with graph's files
        GRAPH_DIR::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "julia"))
        
        # read, prepare and join files with information about edges and their attributes
        graph_structure::DataFrame = read_edges(GRAPH_DIR)
        # read and prepare csv file with nodes attributes
        nodes_attrib::DataFrame = read_nodes_attributes(GRAPH_DIR)

        return graph_structure, nodes_attrib
    end

    function create_graph_structure(graph_structure::DataFrame, nodes_attrib::DataFrame)::graph_frame
        terminal_edges::SubDataFrame = subset(graph_structure, :terminal => x -> x .== true, view=true)
        start_edge::SubDataFrame = subset(graph_structure, :start => x -> x .== true, view=true)
        preterminal_edges::SubDataFrame = subset(graph_structure, :preterminal => x -> x .== true, view=true)
        graph::graph_frame = graph_frame(
            nodes_attrib.ids, # nodes::Vector{Int32}
            Tuple.(Tables.namedtupleiterator(@view nodes_attrib[:, [:x, :y, :z]])), # nodes_coordinates::Vector{Tuple{Float64, Float64, Float64}}
            Tuple.(Tables.namedtupleiterator(@view graph_structure[:, [:source_id, :target_id]])), # edges::Vector{Tuple{Int32, Int32}}
            Tuple.(Tables.namedtupleiterator(terminal_edges[:, [:source_id, :target_id]])), # terminal_edges::Vector{Tuple{Int32, Int32}}
            Tuple.(Tables.namedtupleiterator(start_edge[:, [:source_id, :target_id]])), # start_edge::Tuple{Int32, Int32}
            Tuple.(Tables.namedtupleiterator(preterminal_edges[:, [:source_id, :target_id]])), # preterminal_edges::Vector{Tuple{Int32, Int32}}
            graph_structure.flows, # flows::Vector{Float64}
            graph_structure.volumes, # volumes::Vector{Float64}
            graph_structure.element_ids,
            graph_structure.flow_ids,
            graph_structure.volume_ids
        )
        return graph
    end

    function save_graph!(graph::graph_frame, tree_id::String, graph_id::String)
        JLD2_DIR::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs"))
        jldsave(joinpath(JLD2_DIR, "A.jld2"), graph = graph)
    end

end