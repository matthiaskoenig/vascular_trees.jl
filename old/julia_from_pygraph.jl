module Julia_from_pygraph
"""
Module which contains functions for creating x0 and p vectors specific for each type of PYTHON MERGED graph for ODESolve function.
    These matrices are needed for ODE function from Pharmacokinetic_models.jl.

Also this module can be used for checking ODE functions from Pharmacokinetic_models.jl.
    1. For this uncomment:
        #include("./Pharmacokinetic_models.jl")
        #import .Pharmacokinetic_models: f_dxdt!
        __init__ function
    2. Specialize what graph do you want to use in  __init__ function


Input:
1. pygraph.csv - DataFrame with information about edges of the graph
   and their attributes made in python

Output:
1. x0, p

TODO: adjacency matrix - optimize
"""

# https://juliagraphs.org/Graphs.jl/dev/
using CSV,
    DataFrames,
    Graphs,
    EzXML,
    ParameterizedFunctions,
    GraphDataFrameBridge,
    MetaGraphs,
    OrdinaryDiffEq,
    InteractiveUtils
# using GraphPlot
include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Options: tree_options
# include("./Pharmacokinetic_models.jl")
# import .Pharmacokinetic_models: f_dxdt!

export get_ODE_components

# function __init__()
#     x0, p = get_ODE_components("Rectangle_quad", Int32(10))
#     dx::Vector{Float64} = zeros(size(p, 1))
# end

function get_ODE_components(
    tree_id::String,
    n_node::Int32,
)::Tuple{Vector{Float64},Matrix{Float64}}
    """
    Summarized workflow.
    """
    edges_df, graph = read_graph(tree_id = tree_id, n_node = n_node)
    p = collect_graph_characteristics(edges_df = edges_df, graph = graph)
    x0::Vector{Float64} = zeros(size(p, 1))
    set_initial_values!(x0, p, 1.0)
    return x0, p
end

function read_graph(;
    tree_id::String,
    n_node::Int32,
)::Tuple{DataFrame,MetaDiGraph{Int64,Float64}}
    """
    1. Load and store dataframe with edges information.
    2. Create and store graph from this dataframe.
    3. Return this objects.
    """
    # get graph id for correct path definition
    graph_id::String = "$(tree_id)_$(n_node)"
    # get path to the graph
    GRAPH_PATH::String = normpath(
        joinpath(
            @__FILE__,
            "../../..",
            JULIA_RESULTS_DIR,
            tree_id,
            graph_id,
            "graphs/pygraph.csv",
        ),
    )
    # read csv file with information about edges
    edges_df::DataFrame = DataFrame(CSV.File(GRAPH_PATH))
    # create graph from DataFrame
    graph::MetaDiGraph{Int64,Float64} = MetaDiGraph(
        edges_df,
        :source,
        :target,
        edge_attributes = [:terminal, :start, :group, :is_inflow, :radius, :flow, :length],
    )

    return edges_df, graph
end

function collect_graph_characteristics(;
    edges_df::DataFrame,
    graph::MetaDiGraph{Int64,Float64},
)::Matrix{Float64}
    """
    Workflow for creating correct p vector.
    """

    # calculation of terminal volume
    volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000  # [cm^3] -> [l]
    n_terminals::Int32 = nrow(edges_df[edges_df.terminal, :])
    volume_terminal::Float64 = volume_geometry / n_terminals

    # get full (not sparse) adjacency matrix
    A::Matrix{Float64} = convert(Matrix{Float64}, (adjacency_matrix(graph)))
    # add to the adjacency matrix "edges" connection of terminal node to themselves, 
    # for example, terminal_node_i -> terminal_node_i
    # this is needed for simplifiction of writing ODEs for terminal nodes
    # and for connected to them edges
    modify_adjacency_matrix!(A, graph)

    # get rid of zeros (non existing edges) from adjacency matrix, we do not need them
    edges::Vector{CartesianIndex{2}} = findall(!iszero, A)
    # create parameters matrix - p
    # size of p is determined by: rows - number of connections/edges, columns - number of edges characteristic 
    # that we need to store for writing ODE system correctly
    p::Matrix{Float64} = zeros(length(edges), 8)
    # get and store all the information that we need to write ODE system correctly in the p matrix
    collect_parameters_values!(p, edges, graph, volume_terminal)

    return p
end

function modify_adjacency_matrix!(A::Matrix{Float64}, graph::MetaDiGraph{Int64,Float64})
    @inbounds for (index, _) in pairs(IndexCartesian(), A)
        source_id, target_id = Tuple.(index)
        if (target_id == source_id) && (startswith(props(graph, target_id)[:name], "T_"))
            A[index] = 1.0
        end
    end
end

function collect_parameters_values!(
    p::Matrix{Float64},
    edges::Vector{CartesianIndex{2}},
    graph::MetaDiGraph{Int64,Float64},
    volume_terminal::Float64,
)
    # parameter vector:
    # 1 row = sorce id, traget id, volume value, flow value, whether edge is from inflow system or not,
    # is it source for terminal node, is it target for terminal node, is it a marginal node
    # this order of columns IS CRUCIAL not only here, but also for ODE function
    # create views for convenience
    sources = @view p[:, 1]
    targets = @view p[:, 2]
    volume_values = @view p[:, 3]
    flow_values = @view p[:, 4]
    is_inflow = @view p[:, 5]
    source_for_terminal = @view p[:, 6]
    target_for_terminal = @view p[:, 7]
    is_start_node = @view p[:, 8]
    @inbounds for (ke, edge) ∈ enumerate(edges)
        source_id, target_id = Tuple.(edge)
        sources[ke] = source_id
        targets[ke] = target_id
        # check if an edge is terminal
        if target_id == source_id
            # terminal edge (= terminal node)
            volume_values[ke] = volume_terminal
            is_inflow[ke] = 0.0
            # mark elements connected to terminal nodes
            source_for_terminal[findfirst(x -> x == source_id, targets)] = 1.0
            target_for_terminal[findfirst(x -> x == target_id, sources)] = 1.0
        else
            # other edges
            volume_values[ke] =
                π *
                props(graph, source_id, target_id)[:radius]^2 *
                props(graph, source_id, target_id)[:length] ./ 1000000 # [mm3 --> L]
            flow_values[ke] = props(graph, source_id, target_id)[:flow] # are already converted to [L / min] in python
            is_inflow[ke] = Float64(props(graph, source_id, target_id)[:is_inflow])
            (is_inflow[ke] == 1.0) &&
                (endswith(props(graph, source_id)[:name], "_marginal")) &&
                (is_start_node[ke] = 1.0)
        end
    end
end

function set_initial_values!(
    x0::Vector{Float64},
    p::Matrix{Float64},
    initial_value::Float64,
)
    is_start_node = @view p[:, 8]
    start_nodes_idx = findall(!iszero, is_start_node)
    for start_node_idx in start_nodes_idx
        x0[start_node_idx] = initial_value
    end
end

end
