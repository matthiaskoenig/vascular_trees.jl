module Process_julia_graph
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Definitions: tree_definitions, graph_frame
    import .Utils.Options: graph_options

    include("./processing_helpers.jl")
    import .Processing_helpers: read_edges, read_nodes_attributes, label_special_edges!, create_special_edges!

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
        add_graph_characteristics!(graph_structure)
        graph = create_graph_structure(graph_structure, nodes_attrib, vessel_tree)
        save_graph(graph, tree_id, graph_id, vessel_tree)
    end

#=================================================================================================================================#
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

    function add_graph_characteristics!(graph_structure::DataFrame)
        transform!(graph_structure, [:source_id, :target_id] => ByRow((source_id, target_id) -> ["C_$(source_id)_$(target_id)", "Q_$(source_id)_$(target_id)", "V_$(source_id)_$(target_id)"]) => [:element_ids, :flow_ids, :volume_ids])
        label_special_edges!(graph_structure)
        create_special_edges!(graph_structure)
        assign_ODE_group!(graph_structure)
    end

    function create_graph_structure(graph_structure::DataFrame, nodes_attrib::DataFrame, vessel_tree::String)::graph_frame
        is_inflow::Bool = in(vessel_tree, trees.inflow_trees)
        # reverse edges if the tree is an outflow
        (!is_inflow) && (rename!(graph_structure, [:source_id => :target_id, :target_id => :source_id]))

        terminal_edges = @view graph_structure[graph_structure.terminal .== true, [:source_id, :target_id]]
        start_edge = @view graph_structure[graph_structure.start .== true, [:source_id, :target_id]]
        preterminal_edges = @view graph_structure[graph_structure.preterminal .== true, [:source_id, :target_id]]
        ODE_groups = @view graph_structure[graph_structure.preterminal .== true, :ODE_group]

        # some information must be converted to tuples
        nodes_coordinates::Vector{Tuple{Float64, Float64, Float64}} = Tuple.(Tables.namedtupleiterator(@view nodes_attrib[:, [:x, :y, :z]]))
        edges::Vector{Tuple{Int32, Int32}} = Tuple.(Tables.namedtupleiterator(@view graph_structure[:, [:source_id, :target_id]]))
        terminals::Vector{Tuple{Int32, Int32}} = Tuple.(Tables.namedtupleiterator(terminal_edges))
        start::Vector{Tuple{Int32, Int32}} = Tuple.(Tables.namedtupleiterator(start_edge))
        preterminals::Vector{Tuple{Int32, Int32}} = Tuple.(Tables.namedtupleiterator(preterminal_edges))

        if is_inflow
            pre_elements = zeros(Int64, length(edges))
            post_elements = [[] for _ in eachindex(edges)]
            get_preelements!(pre_elements, edges)
            get_postelements!(post_elements, edges)
        else
            get_outflow_edges_groups!(groups, edges, terminals, preterminals)
        end

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
            ODE_groups,
            pre_elements,
            post_elements
        )
        return graph
    end

    function save_graph(graph::graph_frame, tree_id::String, graph_id::String, vessel_tree::String)
        JLD2_DIR::String = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs"))
        jldsave(joinpath(JLD2_DIR, "$(vessel_tree).jld2"), graph = graph)
    end

#=================================================================================================================================#
    function assign_ODE_group!(graph_structure)
        graph_structure[:, :ODE_group] .= ifelse.(graph_structure[!, :preterminal] .== true, Int16(2),
                                          ifelse.(graph_structure[!, :terminal] .== true, Int16(3),
                                          ifelse.(graph_structure[!, :source_id] .== 0, Int16(0),
                                          Int16(1))))
    end

    function get_preelements!(pre_elements, edges)
        for (ke, element) in enumerate(edges)
            source_id, target_id = element
            for (ked, edge) in enumerate(edges)
                if (edge[2] == source_id != edge[1])  
                    pre_elements[ke] = ked
                    continue
                end
            end
        end
    end

    function get_postelements!(post_elements, edges)
        for (ke, element) in enumerate(edges)
            post_element = []
            source_id, target_id = element
            for (ked, edge) in enumerate(edges)
                if (edge[1] == target_id)
                    push!(post_element, ked)
                end
            end
            append!(post_elements[ke], post_element)
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

end