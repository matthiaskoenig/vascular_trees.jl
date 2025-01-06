module Julia_from_graph
    """
    Module which contains functions for generation f_dxdt function, vector with parameters values 
    and vector with x initials from graph.csv file generated in Python.

    Input:
    1. graph.csv - DataFrame with information about edges of the graph 
       and their attributes

    Output:
    1. x0
    2. p
    3. f_dxdt!
    """
    # https://juliagraphs.org/Graphs.jl/dev/
    using CSV, DataFrames, Graphs, EzXML, ParameterizedFunctions, GraphDataFrameBridge, MetaGraphs, OrdinaryDiffEq, TimerOutputs
    # using GraphPlot
    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Options: graph_options, edge_options

    export get_ODE_components

    function __init__()
        #x0, p, f_dxdt! = get_ODE_components(tree_id="Rectangle_trio", n_node=Int32(10))
        get_ODE_components(tree_id="Rectangle_trio", n_node=Int32(10))
    end

    function get_ODE_components(; tree_id::String, n_node::Int32)
        to = TimerOutput()

        edges_df, graph = read_graph(tree_id=tree_id, n_node=n_node)
        # @timeit to "benchmarking" begin
        A, volume_values, flow_values, group, is_inflow = collect_graph_characteristics(edges_df=edges_df, graph=graph)
        create_f_dxdt(A, volume_values, flow_values, group, is_inflow)
        # end

        # show(to)
        #create_f_dxdt(A=A)

        #elements = collect_edges_metadata(edges_df=edges_df)
        #x0 = collect_initial_values(elements=elements)
        #p = collect_parameters_values(elements=elements)
        
        #f_dxdt! = create_dxdt!(elements=elements, graph=graph)

        #return x0, p, f_dxdt!
    end

    function read_graph(; tree_id::String, n_node::Int32)::Tuple{DataFrame, MetaDiGraph{Int64, Float64}}

        graph_id = "$(tree_id)_$(n_node)"
        GRAPH_PATH = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "graphs/graph.csv"))
        # read csv file with information about edges - main file here
        edges_df = DataFrame(CSV.File(GRAPH_PATH))
        # create graph from DataFrame
        graph = MetaDiGraph(edges_df, :source, :target, edge_attributes=[:terminal, :start, :group, :is_inflow, :radius, :flow, :length])

        return edges_df, graph
    end

    function collect_graph_characteristics(; edges_df::DataFrame, graph::MetaDiGraph{Int64, Float64})::Tuple{Array{Int8}, Array{Float64}, Array{Float64}, Array{String}, Array{Bool}}
        # Return a sparse adjacency matrix for a graph, indexed by [u, v] vertices, so rows - sources, columns - targets
        A::Array{Int8} = Matrix(adjacency_matrix(graph))
        volume_values::Array{Float64} = zeros(size(A))
        flow_values::Array{Float64} = zeros(size(A))
        group = Array{String}(undef, size(A))
        is_inflow = Array{Bool}(undef, size(A))

        volume_geometry::Float64 = (0.100 * 0.100 * 0.10) / 1000  # [cm^3] -> [l]
        n_terminals::Int32 = nrow(edges_df[edges_df.terminal, :])
        volume_terminal::Float64 = volume_geometry / n_terminals

        @inbounds for target_id in 1:size(A, 2), source_id in 1:size(A, 1)
            volume_value = @view volume_values[source_id, target_id]
            if A[source_id, target_id] == 1
                volume_value = π * props(graph, source_id, target_id)[:radius]^2 * props(graph, source_id, target_id)[:length]
                flow_values[source_id, target_id] = props(graph, source_id, target_id)[:flow]
                group[source_id, target_id] = props(graph, source_id, target_id)[:group]
                is_inflow[source_id, target_id] = props(graph, source_id, target_id)[:is_inflow]
            end
            if (target_id == source_id .&& startswith(props(graph, target_id)[:name], "T_"))
                A[source_id, target_id] = 1
                volume_value = volume_terminal
                group[source_id, target_id] = "T"
                is_inflow[source_id, target_id] = false
            end
        end

        return A, volume_values, flow_values, group, is_inflow
    end

    function create_f_dxdt(A, volume_values, flow_values, group, is_inflow)
        #Mustache julia
        #https://juliaci.github.io/PkgTemplates.jl/stable/user/#Custom-Template-Files

        
        elements::Vector{CartesianIndex} = findall(!iszero, A)
        dx::Vector{Float64} = zeros(length(A))
        x::Array{Float64} = zeros(size(A))

        for element ∈ elements
            # retrieve information for element
            source_id, target_id = Tuple.(element)
            element_index = source_id + (target_id - 1) * size(A, 1)

            # species before
            pre_elements = findall(!iszero, @view A[:, source_id]) 

            # species after
            post_elements = findall(!iszero, @view A[target_id, :])

            # marginal inflow element
            if (is_inflow[element] .&& length(pre_elements) == 0) 

                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;

                # species after
                for _ ∈ post_elements
                    dx[element_index] -= flow_values[element] * x[element] / volume_values[element]
                end

            # in -> inflow & inflow element
            elseif (is_inflow[element] .&& length(pre_elements) != 0) 

                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                # or
                # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[element] * x[pre_element, source_id] / volume_values[element]
                end
                
                for post_element ∈ post_elements
                    dx[element_index] -= flow_values[target_id, post_element] * x[element] / volume_values[element]
                end
            
            # terminal element
            elseif (!is_inflow[element] .&& source_id == target_id)

                # dT/dt = QA * A / VT - QV * T / VT
                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[pre_element, source_id] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for post_element ∈ post_elements
                    dx[element_index] -= flow_values[target_id, post_element] * x[element] / volume_values[element]
                end

            # outflow element & outflow -> out
            elseif (!is_inflow[element] .&& source_id != target_id)
                # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
                # or
                # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[pre_element, source_id] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for _ in ∈ post_elements
                    dx[element_index] -= flow_values[element] * x[element] / volume_values[element]
                end

            end

        end

        @show dx
        
    end

    function collect_edges_metadata(; edges_df::DataFrame)::Vector{edge_options}
    
        # gplot(graph)

        element_suffix::String = ""
        element_id::String = ""
        
        flow_value::Float64 = 0.0
        volume_value::Float64 = 0.0

        elements = Vector{edge_options}(nothing, length(nrows(edges_df))+n_terminals)
        initial = 0
        for row in Tables.namedtupleiterator(edges_df)
            initial = 0
            element_suffix = split(row.source, "_")[2] * "_" * split(row.target, "_")[2]
            element_id = row.group * "_" * element_suffix
            flow_value = row.flow
            volume_value = π * row.radius^2 * row.length
            if row.is_inflow & occursin("marginal", element_suffix)
                initial = 10
            end
            push!(elements, edge_options(id=element_id, 
                                         source=row.source, 
                                         target=row.target, 
                                         volume_id="V$(element_id)",
                                         volume_value=volume_value, 
                                         flow_id="Q$(element_id)",
                                         flow_value=flow_value, 
                                         terminal=row.terminal, 
                                         start=row.start, 
                                         inflow=row.is_inflow, 
                                         initial=initial))
            if startswith("T_", row.target)
                push!(elements, edge_options(row.target, "", "", volume_terminal, 0, true, false, false, initial))
            end
        end

        return elements
    end

    function collect_initial_values(; elements::Vector{edge_options})::Vector{Float64}
        return [element.initial for element in elements]
    end

    function collect_parameters_values(; elements::Vector{edge_options})::Vector{Float64}
        volume_values::Vector{Float64} = [element.volume_value for element in elements] 
        flow_values::Vector{Float64} = [element.flow_value for element in elements] 

        return reduce(vcat, (volume_values, flow_values))
    end

    function create_dxdt(; elements::Vector{edge_options}, graph::MetaDiGraph)::Function
        element_index = 0
        pre_element_index = 0
        post_element_index = 0
    
        function f_dxdt(dx::Vector{Float64}, x::Vector{Float64}, p::Vector{Float64}, t::Float64)::Nothing
            for (element_index, element) in enumerate(elements)
                # initialize equation
                dx[element_index]::Float64 = 0.0
                # marginal inflow element
                if element.inflow & occursin("marginal", element.id)
                    # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
                    # species after
                    for _ in outneighbors(graph, graph[:name][element.target])
                        dx[element_index] -= p[length(elements)+element_index] * x[element_index] / p[element_index]
                    end
                # in -> inflow & inflow element
                elseif element.inflow & !occursin("marginal", element.id)
                    # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                    # or
                    # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                    # species before
                    for pre_element in inneighbors(graph, graph[:name][element.source])
                        pre_element_index = findall(x -> x.source == props(graph, pre_element)[:name] .&&  x.target == element.source, elements)[1]
                        dx[element_index] += p[length(elements)+element_index] * x[pre_element_index] / p[element_index]
                    end
                    # species after
                    for post_element in outneighbors(graph, graph[:name][element.target])
                        post_element_index = findall(x -> x.target == props(graph, post_element)[:name] .&&  x.source == element.target, elements)[1]
                        dx[element_index] -= p[length(elements)+post_element_index] * x[element_index] / p[element_index]
                    end
                # outflow element & outflow -> out
                else !element.inflow & !startswith("T_", element.source)
                    # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
                    # or
                    # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
                    # species before
                    for pre_element in inneighbors(graph, graph[:name][element.source])
                        pre_element_index = findall(x -> x.source == props(graph, pre_element)[:name] .&&  x.target == element.source, elements)[1]
                        dx[element_index] += p[length(elements)+pre_element_index] * x[pre_element_index] / p[element_index]
                    end
                    # species after
                    for _ in outneighbors(graph, graph[:name][element.target])
                        dx[element_index] -= p[length(elements)+element_index] * x[element_index] / p[element_index]
                    end
                end
                # inflow -> terminal node
                if startswith("T_", element.target)
                    # dT/dt = QA * A / VT - QV * T / VT
                    # species before
                    for pre_element in inneighbors(graph, graph[:name][element.target])
                        pre_element_index = findall(x -> x.source == props(graph, pre_element)[:name] .&&  x.target == element.target, elements)[1]
                        dx[element_index] += p[length(elements)+pre_element_index] * x[pre_element_index] / p[element_index]
                    end
                    # species after
                    for post_element in outneighbors(graph, graph[:name][element.target])
                        post_element_index = findall(x -> x.target == props(graph, post_element)[:name] .&&  x.source == element.target, elements)[1]
                        dx[element_index] -= p[length(elements)+post_element_index] * x[element_index] / p[element_index]
                    end
                # terminal node -> outflow
                elseif startswith("T_", element.source)
                    # dV/dt = QV * T / VV - QV * V / VV
                    # species before
                    for pre_element in inneighbors(graph, graph[:name][element.target])
                        pre_element_index = findall(x -> x.source == props(graph, pre_element)[:name] .&&  x.target == element.target, elements)[1]
                        dx[element_index] += p[length(elements)+element_index] * x[pre_element_index] / p[element_index]
                    end
                    # species after
                    for _ in outneighbors(graph, graph[:name][element.target])
                        dx[element_index] -= p[length(elements)+element_index] * x[element_index] / p[element_index]
                    end
                end
            end
        end
        return f_dxdt!
    end

end