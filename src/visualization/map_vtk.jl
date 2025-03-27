module Map_VTK
include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Options: graph_options
import .Utils.Definitions: tree_definitions

using DataFrames, Arrow

const trees::tree_definitions = tree_definitions()
# === Graph options ===
# options for graph, i.e., number of nodes and type of tree
const g_options = graph_options(
    n_nodes = [10],  #10
    tree_configurations = [
        "Rectangle_quad",
        # "Rectangle_trio",
    ],
)

coordinates = (0.0, 0.0, 0.0)

# Basic information about the tree that differs between its types (Rectangle_quad, trio, etc.)
# and which is used repeatedly in its processing
Base.@kwdef struct Tree_structure
    tree_configuration::String
    n_node::Integer
    graph_id::String = "$(tree_configuration)_$(n_node)"
    tree_components::Dict{Symbol,Vector{String}} = trees.vascular_trees[tree_configuration]
    vascular_trees::Vector{String} = reduce(vcat, values(tree_components))
    GRAPH_DIR::String = normpath(
        joinpath(@__FILE__, "../../..", JULIA_RESULTS_DIR, tree_configuration, graph_id),
    )
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

for tree_configuration ∈ g_options.tree_configurations
    for n_node in g_options.n_nodes
        tree_info = Tree_structure(; tree_configuration = tree_configuration, n_node = n_node)
        for vascular_tree in ["A"]
            coord_start::Bool = false
            GRAPH_PATH::String = joinpath(tree_info.GRAPH_DIR, "graphs", "$(vascular_tree).arrow")
            VTK_PATH::String = joinpath(tree_info.GRAPH_DIR, "julia", "$(vascular_tree).vtk")

            p = get_graph_parameters(GRAPH_PATH)

            is_inflow = p[2]
            nodes_ids = p[4]
            species_ids = p[6]
            flows_values = p[7]

            nodes_coordinates = [round.(x; digits=8) for x in p[5]]
            species_to_nodes_ids = Vector{String}(undef, count(!ismissing, nodes_coordinates))
            flows_to_node_ids = Vector{AbstractFloat}(undef, count(!ismissing, nodes_coordinates))

            nodes_lines = Vector{Integer}(undef, count(!ismissing, nodes_coordinates))
            
            
            if is_inflow
                equality_edge_to_node = 2
            else
                equality_edge_to_node = 1
            end

            open(VTK_PATH) do f
                # line_number
                line = 0
                # read till end of file
                while ! eof(f)
                    # read a new / next line for every iteration           
                    content = readline(f)
                    if coord_start && isempty(content)
                        coord_start = false
                    end
                    if coord_start
                        coordinates = Tuple(parse.(Float64, split(content, " ")))
                        for (kn, node_coordinates) in enumerate(skipmissing(nodes_coordinates))
                            if isequal(node_coordinates, coordinates)
                                nodes_lines[kn] = line
                                for (ke, edge_id) in enumerate(p[3])
                                    if edge_id[1] != edge_id[2]
                                        if edge_id[equality_edge_to_node] == nodes_ids[kn]
                                            species_to_nodes_ids[kn] = species_ids[ke]
                                            flows_to_node_ids[kn] = flows_values[ke]
                                            continue
                                        end
                                    end
                                end
                                continue
                            end
                        end     
                    end
                    if startswith(content, "POINTS")
                        coord_start = true    
                    end      
                    line += 1
                end
            end
            map_df = DataFrame(
                node_line = nodes_lines,
                node_id = skipmissing(nodes_ids),
                species_id = species_to_nodes_ids,
                flows_values = flows_to_node_ids
            )
            @show map_df
        end
    end
end



#duplicate function from julia_from_jgraph.jl
function get_graph_parameters(GRAPH_PATH::String, n_inflow::Integer)
    graph = DataFrame(Arrow.Table(GRAPH_PATH))
    p = terminal_parameters(; id = "T", x_affiliations = Array{String}(graph[1:(n_inflow+1), :]), flow_values = Array{Float64}(graph[(n_inflow+2):(n_inflow+2+n_inflow), :]) .* flow_scaling_factor, volumes = Float64((graph[end, 1])))

    return p
end


end