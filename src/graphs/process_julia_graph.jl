module Process_Julia_Graph
"""
Module for processing julia graph files from SyntheticVascularTrees.jl.
    Only works with "Rectangle_quad" and "Rectangle_trio". Otherwise "tree_definitions" must be modified in Utils.jl.

Input: .lg and .csv files for every individual vessel tree (arterial, portal, etc.)

Output: .arrow (table) for every individual vessel tree (arterial, portal, etc.) + 
    one table for terminal nodes which serve as connection points between the trees.

Idea of this module: to get, to store, and to save all the information that we need for correct ODEs from the graph. 

TODO: Optimize code
"""

include("../utils.jl")
import .Utils.Options: graph_options
import .Utils.Definitions: tree_definitions
import .Utils: JULIA_RESULTS_DIR
include("./processing_helpers.jl")

include("./process_individual_trees.jl")
import .Process_Individual_Trees: process_julia_graph
include("./process_terminal_nodes.jl")
import .Process_Terminal_Nodes: process_terminal_nodes

# Already specified in utils.jl
const trees::tree_definitions = tree_definitions()

# === Graph options ===
# options for groph, i.e., number of nodes and type of tree
g_options = graph_options(
    n_nodes = [100],  #750, 1000, 1250, 1500
    tree_configurations = ["Rectangle_quad"],
)

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

for tree_configuration ∈ g_options.tree_configurations
    for n_node ∈ g_options.n_nodes
        tree_info =
            Tree_structure(; tree_configuration = tree_configuration, n_node = n_node)
        process_julia_graph(tree_info)
        process_terminal_nodes(tree_info)
    end
end

end
