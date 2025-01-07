module Julia_from_graph
    """
    Module which contains functions for generation f_dxdt function, vector with parameters values 
    and vector with x initials from graph.csv file generated in Python.

    Input:
    1. graph.csv - DataFrame with information about edges of the graph 
       and their attributes

    Output:
    1. f_dxdt
    """
    # https://juliagraphs.org/Graphs.jl/dev/
    using CSV, DataFrames, Graphs, EzXML, ParameterizedFunctions, GraphDataFrameBridge, MetaGraphs, OrdinaryDiffEq, TimerOutputs
    # using GraphPlot
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Options: graph_options, edge_options

    export get_ODE_components

    function __init__()
        get_ODE_components(tree_id="Rectangle_trio", n_node=Int32(10))
    end

    function get_ODE_components(; tree_id::String, n_node::Int32)
        to = TimerOutput()

        edges_df, graph = read_graph(tree_id=tree_id, n_node=n_node)
        @timeit to "benchmarking" begin
            p = collect_graph_characteristics(edges_df=edges_df, graph=graph)
            A = @view p[:, 1:Int8(size(p, 2)/4)]
            dx::Vector{Float64} = zeros(length(A))
            x::Array{Float64} = zeros(size(A))
            f_dxdt!(dx, x, p)
        end

        show(to)
        #create_f_dxdt(A=A)

        #elements = collect_edges_metadata(edges_df=edges_df)
        #x0 = collect_initial_values(elements=elements)
        #p = collect_parameters_values(elements=elements)
        
        #f_dxdt! = create_dxdt!(elements=elements, graph=graph)

        #return x0, p, f_dxdt!
    end

    function f_dxdt!(dx::Vector{Float64}, x::Array{Float64}, p::Array{Float64})

        # create views for convenience
        # p contains adjacency matrix, volumes, flows, 
        # is the edge from inflow system or not
        A = @view p[:, 1:Int8(size(p, 2)/4)]
        volume_values = @view p[:, size(A, 2)+1:size(A, 2)*2]
        flow_values = @view p[:, size(A, 2)*2+1:size(A, 2)*3]
        is_inflow = @view p[:, size(A, 2)*3+1:size(p, 2)]

        # take only existing edges from the adjacency matrix
        elements::Vector{CartesianIndex} = findall(!iszero, A)

        for element ∈ elements
            # retrieve information for element
            source_id, target_id = Tuple.(element)
            element_index = source_id + (target_id - 1) * size(A, 1)
            # species before the slement
            pre_elements = findall(!iszero, @view A[:, source_id]) 
            # species after the slement
            post_elements = findall(!iszero, @view A[target_id, :])

            # distinguish elements connected to terminal nodes
            terminal_source::Bool = false
            if (length(post_elements) == 1) && (target_id == post_elements[1])
                terminal_source = true
            end
            terminal_target::Bool = false
            if (length(pre_elements) == 1) && (source_id == pre_elements[1])
                terminal_target = true
            end

            # write equations
            # marginal inflow element & in -> inflow & inflow element not connected to terminal
            if (terminal_source == false .&& is_inflow[element] == 1.0) 
                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
                # or
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                # or
                # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[element] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for post_element ∈ post_elements
                    dx[element_index] -= flow_values[target_id, post_element] * x[element] / volume_values[element]
                end

            # inflow element connected to terminal
            elseif (terminal_source == true .&& is_inflow[element] == 1.0 .&& length(pre_elements) != 0) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[element] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for _ ∈ post_elements
                    dx[element_index] -= flow_values[element] * x[element] / volume_values[element]
                end
            
            # terminal element
            elseif (is_inflow[element] == 0.0 .&& source_id == target_id)
                # dT/dt = QA * A / VT - QV * T / VT

                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[pre_element, source_id] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for post_element ∈ post_elements
                    dx[element_index] -= flow_values[target_id, post_element] * x[element] / volume_values[element]
                end

            # outflow element connected to terminal
            elseif (terminal_target == true .&& is_inflow[element] == 0.0 .&& source_id != target_id)
                # dV/dt = QV * T / VT - QV * V / VV;

                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[element] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for _ ∈ post_elements
                    dx[element_index] -= flow_values[element] * x[element] / volume_values[element]
                end

            # outflow element & outflow -> out not connected to terminal
            elseif (terminal_target == false .&& is_inflow[element] == 0.0 .&& source_id != target_id)
                # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
                # or
                # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[pre_element, source_id] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for _ ∈ post_elements
                    dx[element_index] -= flow_values[element] * x[element] / volume_values[element]
                end

            end

        end

        #@show dx
        
    end


    function read_graph(; tree_id::String, n_node::Int32)::Tuple{DataFrame, MetaDiGraph{Int64, Float64}}

        # get graph id
        graph_id = "$(tree_id)_$(n_node)"
        # get path to the graph
        GRAPH_PATH = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs/graph.csv"))
        # read csv file with information about edges
        edges_df = DataFrame(CSV.File(GRAPH_PATH))
        # create graph from DataFrame
        graph = MetaDiGraph(edges_df, :source, :target, 
                            edge_attributes=[:terminal, :start, :group, :is_inflow, :radius, :flow, :length])

        return edges_df, graph
    end

    function collect_graph_characteristics(; edges_df::DataFrame, graph::MetaDiGraph{Int64, Float64})::Array{Float64}
        
        # create a parameter vector, which will contain:
        # adjacency matrix, volume values, flow_values, whether edge is from inflow system or not
        # THIS ORDER IS CRUCIAL
        p = repeat(convert(Matrix{Float64}, (adjacency_matrix(graph))), outer=(1, 4))

        # create views for convenience
        # adjacency matrix for a graph, indexed by [u, v] vertices, so rows - sources, columns - targets
        A = @view p[:, 1:Int8(size(p, 2)/4)]
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


        for (index, value) in pairs(IndexCartesian(), A)
            source_id, target_id = Tuple.(index)
            if value == 1.0
                volume_values[index] = π * props(graph, source_id, target_id)[:radius]^2 * props(graph, source_id, target_id)[:length]
                flow_values[index] = props(graph, source_id, target_id)[:flow]
                is_inflow[index] = Float64(props(graph, source_id, target_id)[:is_inflow])
            end
            if (target_id == source_id .&& startswith(props(graph, target_id)[:name], "T_"))
                value = 1.0
                volume_values[index] = volume_terminal
                is_inflow[index] = 0.0
            end
        end

        #@show p

        return p
    end

end