module Map_VTK

"""


Note for understanding the workflow with T.arrow file (terminal part of the tree): 
    each tree component has terminal nodes. So, when we have terminal node 1 (T_1),
    this means that we have for it instance (affiliation) from arterial tree, venous tree, etc.
    Each of this instance has DIFFERENT coordinates, so we need to store them all.
    Example for Rectangle quad: one terminal node from T.arrow has four coordinate points.
"""



include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Options: tree_options
import .Utils.Definitions: tree_definitions

using DataFrames, Arrow, CSV

const trees::tree_definitions = tree_definitions()
# === Graph options ===
# options for graph, i.e., number of nodes and type of tree
const t_options = tree_options(
    n_nodes = [10, 100, 1000],  #10
    tree_configurations = [
        "Rectangle_quad",
        # "Rectangle_trio",
    ],
)

# Basic information about the tree that differs between its types (Rectangle_quad, trio, etc.)
# and which is used repeatedly in its processing
Base.@kwdef struct Tree_structure
    tree_configuration::String
    n_node::Integer
    graph_id::String = "$(tree_configuration)_PVL_ligated_$(n_node)"
    tree_components::Dict{Symbol,Vector{String}} = trees.vascular_trees[tree_configuration]
    vascular_trees::Vector{String} = reduce(vcat, values(tree_components))
    GRAPH_DIR::String = normpath(
        joinpath(@__FILE__, "../../..", JULIA_RESULTS_DIR, tree_configuration, graph_id),
    )
end


function __init__()
    for tree_configuration ∈ t_options.tree_configurations
        for n_node in t_options.n_nodes
            tree_info = Tree_structure(; tree_configuration = tree_configuration, n_node = n_node)
            map_vtk_with_jgraphs(tree_info)
        end
    end
end    

function map_vtk_with_jgraphs(tree_info::Tree_structure)
    """
    Full workflow of mapping one full tree with its vtk files.
    """
    # load arrow file with terminal part info
    # make first mappings for this part
    terminal_nodes_mapping = map_vtk_with_terminal_file(tree_info)
    for vascular_tree in tree_info.vascular_trees
        # map corresponding arrow files with individual tree components arrow files
        # + finish mapping for terminal part 
        map_vtk_with_individual_jgraphs!(terminal_nodes_mapping, tree_info, vascular_tree)
    end
    # create dataframe with mapping info for TERMINAL part
    df_map_terminal_file = DataFrame(
        :file_name => terminal_nodes_mapping.file_names,
        :node_line => terminal_nodes_mapping.node_line,
        :species_ids => terminal_nodes_mapping.species_ids)
    # Write DataFrame to CSV
    CSV.write(joinpath(tree_info.GRAPH_DIR, "graphs", "T_map_vtk.csv"), df_map_terminal_file)
end

function map_vtk_with_terminal_file(tree_info::Tree_structure)::NamedTuple
    """
    Function that loads T.arrow file, gets parameters that we need for mapping with vtk from this file,
        make first mapping:
        1. Terminal node id from this file to the corresponding vascular trees 
            (all needed info for this is already in T.arrow file)
        + gets terminal nodes coordinates, so we can map further terminal nodes from T.arrow to
        corresponding nodes in each vascular tree arrow file based on the same coordinates.
        + prepare vector which will further contain nodes lines from corresponding vascular tree arrow files.
    """
    # get path to the T.arrow file
    TERMINAL_NODES_PATH::String = joinpath(tree_info.GRAPH_DIR, "graphs", "T.arrow")

    # get number of inflows and number of tree components (number of vascular systems)
    # they are needed to separate T.arrow file on its components correctly
    n_inflows::Integer = length(tree_info.tree_components[:inflow_trees])
    n_tree_components::Integer = length(tree_info.tree_components[:inflow_trees]) + length(tree_info.tree_components[:outflow_trees])

    # get info needed for mapping
    terminal_info = get_graph_parameters(TERMINAL_NODES_PATH, n_inflows, n_tree_components)


    # =========== Get more needed further information
    # number of terminal nodes
    n_terminals = size(terminal_info.x_affiliations)[2]
    # ids of terminal species
    terminal_ids = view(terminal_info.x_affiliations, 1, :)
    # get ids of NON TERMINAL species from T.arrow file. They are all from inflow
    # tree components
    inflow_species_ids = vec(view(terminal_info.x_affiliations, 2:size(terminal_info.x_affiliations)[1], :))
    # get ids of inflow tree components
    inflow_file_names = map(x-> split(x, "_")[1], inflow_species_ids)
    
    # allocate vectors needed below
    species_ids = Vector{String}(undef, n_terminals * (n_tree_components + n_inflows))
    file_names = Vector{String}(undef, n_terminals * (n_tree_components + n_inflows))
    node_line = Vector{Integer}(undef, n_terminals * (n_tree_components + n_inflows))
    nodes_coordinates = []

    # terminal node's index, we start from the first one
    terminal_node_idx = 1
    # ODE for each terminal nodes includes several species, these are their index
    terminal_component_idx = 1
    # iterate over terminal coordinates
    for (kc, terminal_node_coordinates) in enumerate(terminal_info.terminal_coordinates)
        # store terminal ids in a common for all species from T.arrow file vector
        species_ids[kc] = terminal_ids[terminal_node_idx]
        # store terminal coordinates in a common for all coordinates from T.arrow file vector
        push!(nodes_coordinates, terminal_node_coordinates)
        # store corresponding for these terminal_node_coordinates tree component id
        # tree component id is equal to the vtk file name
        file_names[kc] = terminal_info.terminal_coordinates_affiliations[kc]

        # update iterators
        terminal_component_idx += 1
        if terminal_component_idx == n_tree_components+1
            terminal_component_idx = 1
            terminal_node_idx += 1
        end
    end
    # store ids of non terminal species in a common for all species from T.arrow file vector
    species_ids[n_terminals*n_tree_components+1:end] = inflow_species_ids
    # store tree component ids for non terminal species in a common vector
    # so we can say from which tree component this non terminal species is
    file_names[n_terminals*n_tree_components+1:end] = inflow_file_names

    # store all colelcted info in tuple
    terminal_nodes_mapping = (
        species_ids = species_ids,
        file_names = file_names,
        terminal_node_coordinates = [round.(node_coordinates; digits=8) for node_coordinates in nodes_coordinates],
        node_line = node_line
    )

    return terminal_nodes_mapping
end

function map_vtk_with_individual_jgraphs!(terminal_nodes_mapping, tree_info, vascular_tree)
    # get paths to arrow, vtk files of the tree component (vascular tree) and get path to where the csv file with mapping will be stored
    GRAPH_PATH, VTK_PATH, MAPPING_RESULTS_PATH = paths_initialization(tree_info.GRAPH_DIR, vascular_tree)
    # load and store info needed for mapping from corresponding vascular tree arrow file
    graph_info = get_graph_parameters(GRAPH_PATH)

    # create these variable for convenience
    is_inflow = graph_info[2]
    nodes_ids = graph_info[4]
    species_ids = graph_info[6]
    flows_values = graph_info[7]

    # nodes coordinates from arrow file should be rounded,
    # so we could compare them to the coordinates from vtk file
    nodes_coordinates = [round.(x; digits=8) for x in graph_info[5]]

    # allocate some needed below vectors
    # we MUST be sure that components of these vectors with the same ids BELONG TO THE SAME NODE
    species_to_nodes_ids = Vector{String}(undef, count(!ismissing, nodes_coordinates))
    flows_to_node_ids = Vector{AbstractFloat}(undef, count(!ismissing, nodes_coordinates))
    nodes_lines = Vector{Integer}(undef, count(!ismissing, nodes_coordinates))
    
    # species ids in the model belong to edges, for the visualization we need to map species 
    # to nodes from vtk file 
    # for the inflow this mapping is done by TARGET NODE ID from the edge (index 2 in the edge tuple) 
    # must be equal to NODE ID
    # for the outflow: SOURCE NODE ID (index 1 in the edge tuple) from the edge must be equal to NODE ID
    # create corresponding variable for this
    if is_inflow
        equality_edge_to_node = 2
    else
        equality_edge_to_node = 1
    end

    # variable for understanding whether we are reading paert of the file
    # with coordinates data or not
    is_coord_part::Bool = false
    # open vtk file and read it line by line
    open(VTK_PATH) do vtk_file
        # line_number
        line = 0
        # read till end of file
        while ! eof(vtk_file)
            # read a new / next line for every iteration           
            content = readline(vtk_file)
            # checking whether part with coordinates has ended
            if is_coord_part && isempty(content)
                is_coord_part = false
            end
            # =================== Mapping INDIVIDUAL trees
            # if we are in coordinates part do:
            if is_coord_part
                # split string line into three parts and convert string to float
                # -> this is one point with x, y, z coordinates
                coordinates = Tuple(parse.(Float64, split(content, " ")))
                # find the node with the same coordinates in the info from arrow file, for this:
                for (kn, node_coordinates) in enumerate(skipmissing(nodes_coordinates))
                    # iterate over coordinates from arrow file
                    # kn - index of the node coordinates from arrow file,
                    # all the mapping info for this node will be stored in corresponding vectors
                    # at the index kn
                    # find the same coordinates from vtk and arrow files
                    # same coordinates mean the same node
                    if isequal(node_coordinates, coordinates)
                        # store the line number with the coordinates if this node in vtk file
                        # so we map a node from arrow file to its line in vtk file
                        nodes_lines[kn] = line
                        # we also need to map nodes with species ids, for this:
                        # iterate over the edges of this tree component (because species ids beling to edges)
                        for (ke, edge_id) in enumerate(graph_info[3])
                            # we do not need terminal edges, check for them
                            if edge_id[1] != edge_id[2]
                                # here we check whether the edge can be mapped with this node
                                # the rule for this is written above
                                if edge_id[equality_edge_to_node] == nodes_ids[kn]
                                    # map the edge with the node
                                    species_to_nodes_ids[kn] = species_ids[ke]
                                    # get flow value for this edge/node
                                    flows_to_node_ids[kn] = flows_values[ke]
                                    continue
                                end
                            end
                        end
                        continue
                    end
                end
                # ============== Continuing mapping the TERMINAL part
                # terminal part is artificial, it consists of the parts of every individual tree
                # so terminal part does not have a separate vtk file, but
                # we need to map nodes from terminal file to the vtk files of individual trees
                # here we map only TERMINAL NODES, other nodes are mapped below
                for (kt, terminal_nodes_coordinates) in enumerate(terminal_nodes_mapping.terminal_node_coordinates)
                    # mapping is based on: same coordinates in the vtk file and in file for terminal part
                    # and these coordinates should come from the SAME tree component
                    if isequal(terminal_nodes_coordinates, coordinates) .&& (terminal_nodes_mapping.file_names[kt] == vascular_tree)
                        # store the line nu,ber from vtk
                        terminal_nodes_mapping.node_line[kt] = line
                        continue
                    end
                end
            end
            # check whether the part of vtk with coordinates has started
            if startswith(content, "POINTS")
                is_coord_part = true    
            end      
            # update line number
            line += 1
        end
    end
    # ============== Continuing mapping the TERMINAL part
    # here we map NON TERMINAL species from terminal file to vtk 
    # we already has done it with species from INDIVIDUAL TREES (mapped them with vtk),
    # so we can find the same species from terminal file and from individual tree file
    # and get the corresponding line number from vtk
    for (ks, species_id) in enumerate(species_to_nodes_ids)
        for (kst, species_id_from_terminal_file) in enumerate(terminal_nodes_mapping.species_ids)
            # if the species id from terminal file is the same to the species id from the individual vascular tree file
            # store he mapping line number from vtk
            if species_id == species_id_from_terminal_file
                terminal_nodes_mapping.node_line[kst] = nodes_lines[ks]
                continue
            end
        end
    end
    # create dataframe with mapping info for INDIVIDUAL tree component
    map_df = DataFrame(
        file_id = [vascular_tree for _ in 1:length(nodes_lines)],
        node_line = nodes_lines,
        node_id = [node_id for node_id in skipmissing(nodes_ids)],
        species_id = species_to_nodes_ids,
        flows_values = flows_to_node_ids,
        nodes_coordinates = [node_coordinates for node_coordinates in skipmissing(nodes_coordinates)]
    )
    # Write DataFrame to CSV
    CSV.write(MAPPING_RESULTS_PATH, map_df)
end


function paths_initialization(GRAPH_DIR::String, vascular_tree::String)::Tuple{String, String, String}
    """
    Function that returns paths to arrow file of this vascular tree, to its vtk file
        and to where the mapping results will be stored.  
    """
    GRAPH_PATH::String = joinpath(GRAPH_DIR, "graphs", "$(vascular_tree).arrow")
    VTK_PATH::String = joinpath(GRAPH_DIR, "julia", "$(vascular_tree).vtk")
    MAPPING_RESULTS_PATH::String = joinpath(GRAPH_DIR, "graphs", "$(vascular_tree)_map_vtk.csv")
    return GRAPH_PATH, VTK_PATH, MAPPING_RESULTS_PATH
end


#duplicate function from julia_from_jgraph.jl
function get_graph_parameters(GRAPH_PATH::String)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    p = (
        graph.vascular_tree_id[1],
        graph.is_inflow[1],
        Vector(graph.all_edges),
        Vector(graph.nodes_ids),
        Vector(graph.nodes_coordinates),
        Vector(graph.species_ids),
        Vector(graph.flows),
    )
    return p
end


#duplicate function from julia_from_jgraph.jl
function get_graph_parameters(GRAPH_PATH::String, n_inflow::Integer, n_tree_components::Integer)
    """"
    Note for understanding (Duplicate): each tree component has terminal nodes. So, when we have terminal node 1 (T_1),
        this means that we have for it instance (affiliation) from arterial tree, venous tree, etc.
        Each of this instance has DIFFERENT coordinates, so we need to store them all.
        Example for Rectangle quad: one terminal node from T.arrow has four coordinate points.
    :id - terminal part id
    :x_affiliations - species ids
    :terminal_coordinates - coordinates of only TERMINAL NODES (only for Ts, not for preterminals)
    :terminal_coordinates_affiliations - array for understanding what coordinate point came from what tree component.
    """
    # load arrow file
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    # get info that we need for mapping for terminal file
    p = (id = "T", 
        x_affiliations = Array{String}(graph[1:(n_inflow+1), :]), 
        terminal_coordinates = Array(graph[(3*n_inflow+4):(3*n_inflow+3+n_tree_components), :]), 
        terminal_coordinates_affiliations = Array(graph[(3*n_inflow+4+n_tree_components):(3*n_inflow+3+2*n_tree_components), :]))
    return p
end


end