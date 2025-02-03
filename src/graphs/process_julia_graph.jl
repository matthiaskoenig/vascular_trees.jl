module Process_julia_graph
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Definitions: tree_definitions, graph_frame
    import .Utils.Options: graph_options

    include("./processing_helpers.jl")
    import .Processing_helpers: read_edges, read_nodes_attributes

    using DataFrames, JLD2, InteractiveUtils

    """
    Module for processing julia graph files from SyntheticVascularTrees.jl.
        Only works with "Rectangle_quad" and "Rectangle_trio". Otherwise "tree_definitions" must be modified in Utils.jl.

    Input: .lg and .csv files for every individual vessel tree (arterial, portal, etc.)

    Output: .JLD2 with graph_frame for every individual vessel tree (arterial, portal, etc.)

    Idea of this module: to get, to store, and to save all the information that we need for correct ODEs from the graph. 

    TODO: Get rid of zeros in DataFrame
    """

    # ============ Specify options
    g_options::graph_options = graph_options(
        n_nodes=[10],  #10, 30, 50, 100, 250, 500, 750, 1000, 1250, 1500, 1750   300000, 400000, 500000, 1000000
        tree_ids=[
            "Rectangle_quad",
            # "Rectangle_trio",
            ],
    )

    # Already specified in utils.jl
    trees::tree_definitions = tree_definitions()

    function __init__()
        for tree_id ∈ g_options.tree_ids, n_node ∈ g_options.n_nodes
            process_julia_graph(tree_id, n_node)
        end
    end

    # Main function: workflow for whole tree
    function process_julia_graph(tree_id::String, n_node::Int32)
        # get graph id for correct path definition
        graph_id::String = "$(tree_id)_$(n_node)"
        println()
        printstyled("------------------------------------------------------------------------------------\n"; color = 124)
        printstyled("   Processing $(graph_id)_$(tree_id)   \n"; color = 9)
        printstyled("------------------------------------------------------------------------------------\n"; color = 124)
        # get directory with graph's files
        GRAPH_DIR::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "julia"))
        for vessel_tree ∈ trees.vascular_trees[tree_id]
            process_individual_tree(GRAPH_DIR, vessel_tree, tree_id, graph_id)
        end
    end

    # Workflow for individual vessel trees
    function process_individual_tree(GRAPH_DIR::String, vessel_tree::String, tree_id::String, graph_id::String)
        GRAPH_PATH, EDGES_PATH, NODES_PATH = paths_initialization(GRAPH_DIR, vessel_tree)
        graph_structure, nodes_attrib = read_graph(GRAPH_PATH, EDGES_PATH, NODES_PATH)
        graph = create_graph_structure(graph_structure, nodes_attrib, vessel_tree)
        save_graph!(graph, tree_id, graph_id, vessel_tree)
    end

    function paths_initialization(GRAPH_DIR::String, vessel_tree::String)::Tuple{String, String, String}
        GRAPH_PATH::String = joinpath(GRAPH_DIR, "$(vessel_tree).lg")
        EDGES_PATH::String = joinpath(GRAPH_DIR, "$(vessel_tree)_edges.csv")
        NODES_PATH::String = joinpath(GRAPH_DIR, "$(vessel_tree)_nodes.csv")
        return GRAPH_PATH, EDGES_PATH, NODES_PATH
    end

    function read_graph(GRAPH_PATH::String, EDGES_PATH::String, NODES_PATH::String)::Tuple{DataFrame, DataFrame}
        # read, prepare and join files with information about edges and their attributes
        graph_structure::DataFrame = read_edges(GRAPH_PATH, EDGES_PATH)
        # read and prepare csv file with nodes attributes
        nodes_attrib::DataFrame = read_nodes_attributes(NODES_PATH)

        return graph_structure, nodes_attrib
    end

    function create_graph_structure(graph_structure::DataFrame, nodes_attrib::DataFrame, vessel_tree::String)::graph_frame
        is_inflow::Bool = in(vessel_tree, trees.inflow_trees)
        # reverse edges if the tree is an outflow
        (!is_inflow) && (rename!(graph_structure, [:source_id => :target_id, :target_id => :source_id]))

        terminal_edges::SubDataFrame = subset(graph_structure, :terminal => x -> x .== true, view=true)
        start_edge::SubDataFrame = subset(graph_structure, :start => x -> x .== true, view=true)
        preterminal_edges::SubDataFrame = subset(graph_structure, :preterminal => x -> x .== true, view=true)

        # some information must be converted to tuples
        nodes_coordinates::Vector{Tuple{Float64, Float64, Float64}} = Tuple.(Tables.namedtupleiterator(@view nodes_attrib[:, [:x, :y, :z]]))
        edges::Vector{Tuple{Int32, Int32}} = Tuple.(Tables.namedtupleiterator(@view graph_structure[:, [:source_id, :target_id]]))
        terminals::Vector{Tuple{Int32, Int32}} = Tuple.(Tables.namedtupleiterator(terminal_edges[:, [:source_id, :target_id]]))
        start::Vector{Tuple{Int32, Int32}} = Tuple.(Tables.namedtupleiterator(start_edge[:, [:source_id, :target_id]]))
        preterminals::Vector{Tuple{Int32, Int32}} = Tuple.(Tables.namedtupleiterator(preterminal_edges[:, [:source_id, :target_id]]))
        groups = Vector{Int16}(undef, length(edges))

        is_inflow ? get_inflow_edges_groups!(groups, edges, terminals, preterminals) : get_outflow_edges_groups!(groups, edges, terminals, preterminals)

        # collect all needed information and store it in one structure
        graph::graph_frame = graph_frame(
            vessel_tree, # vascular_tree_id
            is_inflow, # is_inflow
            nodes_attrib.ids, # nodes::Vector{Int32}
            nodes_coordinates,
            edges, # 
            terminals, # terminal edges
            start, # start edge
            preterminals, # preterminal edges
            graph_structure.flows, # flows::Vector{Float64}
            graph_structure.volumes, # volumes::Vector{Float64}
            graph_structure.element_ids,
            graph_structure.flow_ids,
            graph_structure.volume_ids,
            groups 
        )
        return graph
    end

    function get_inflow_edges_groups!(groups::Vector{Int16}, edges::Vector{Tuple{Int32, Int32}}, terminals::Vector{Tuple{Int32, Int32}}, preterminals::Vector{Tuple{Int32, Int32}})
        @inbounds for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            is_preterminal = in(element, preterminals)
            is_terminal = in(element, terminals)
            # in -> inflow & inflow element not connected to terminal
            if (!is_preterminal .&& !is_terminal .&& source_id != 0)
                groups[ke] = Int16(1)
            # inflow element connected to terminal
            elseif (is_preterminal)
                groups[ke] = Int16(2)
            # terminal element
            elseif (is_terminal)
                groups[ke] = Int16(3)
            # marginal (input) element
            elseif source_id == 0
                groups[ke] = Int16(0)
            end
        end
    end

    function get_outflow_edges_groups!(groups::Vector{Int16}, edges::Vector{Tuple{Int32, Int32}}, terminals::Vector{Tuple{Int32, Int32}}, preterminals::Vector{Tuple{Int32, Int32}})
        @inbounds for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            is_preterminal = in(element, preterminals)
            is_terminal = in(element, terminals)

            # outflow -> out & outflow element not connected to terminal
            if (!is_preterminal .&& !is_terminal .&& target_id != 0)
                groups[ke] = Int16(1)
            # outflow element connected to terminal
            elseif (is_preterminal) 
                groups[ke] = Int16(2)
            # terminal element
            elseif (is_terminal)
                groups[ke] = Int16(3)
            # marginal (outflow) element
            elseif target_id == 0
                groups[ke] = 0
            end
        end
    end

    function save_graph!(graph::graph_frame, tree_id::String, graph_id::String, vessel_tree::String)
        JLD2_DIR::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs"))
        jldsave(joinpath(JLD2_DIR, "$(vessel_tree).jld2"), graph = graph)
    end

end