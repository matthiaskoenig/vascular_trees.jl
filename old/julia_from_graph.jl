# https://juliagraphs.org/Graphs.jl/dev/
using CSV,
    DataFrames,
    Graphs,
    EzXML,
    ParameterizedFunctions,
    GraphDataFrameBridge,
    MetaGraphs,
    DifferentialEquations
using GraphPlot

# read csv file with information about edges - main file here
edges_df = DataFrame(
    CSV.File(
        "/home/mariia/liver_vascularisation/liver_vascular_tree/results/julia_vessel_trees/Rectangle_quad/Rectangle_quad_10/graphs/graph.csv",
    ),
)
# create graph fro# outflow -> outm DataFrame
graph = MetaDiGraph(
    edges_df,
    :source,
    :target,
    edge_attributes = [:terminal, :start, :group, :is_inflow],
)
# gplot(graph)
volume_geometry = (0.100 * 0.100 * 0.10) / 1000  # [cm^3] -> [l]
n_terminals = nrow(edges_df[edges_df.terminal, :])
volume_terminal = volume_geometry / n_terminals


element_suffix = ""
element_id = ""
element_index = 0
pre_element_index = 0
pre_flow_value = 0.0
post_element_index = 0
post_flow_value = 0.0
flow_value = 0.0
volume_value = 0.0

struct Edge
    id::String
    source::String
    target::String
    volume_value::Float32
    flow_value::Float32
    terminal::Bool
    start::Bool
    inflow::Bool
    initial::Float32
end

elements = []
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
    push!(
        elements,
        Edge(
            element_id,
            row.source,
            row.target,
            volume_value,
            flow_value,
            row.terminal,
            row.start,
            row.is_inflow,
            initial,
        ),
    )
    if startswith("T_", row.target)
        push!(
            elements,
            Edge(row.target, "", "", volume_terminal, 0, true, false, false, initial),
        )
    end
end

function vessel_tree!(du, u, p, t)
    for (element_index, element) in enumerate(elements)
        # initialize equation
        du[element_index] = 0
        # marginal inflow element
        if element.inflow & occursin("marginal", element.id)
            # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
            # species after
            for _ in outneighbors(graph, graph[:name][element.target])
                du[element_index] -=
                    element.flow_value * u[element_index] / element.volume_value
            end
            # in -> inflow & inflow element
        elseif element.inflow & !occursin("marginal", element.id)
            # dA/dt = QA * A_pre / VA - Q_post * A / VA;
            # or
            # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
            # species before
            for pre_element in inneighbors(graph, graph[:name][element.source])
                pre_element_index = findall(
                    x ->
                        x.source == props(graph, pre_element)[:name] .&&
                        x.target == element.source,
                    elements,
                )[1]
                du[element_index] +=
                    element.flow_value * u[pre_element_index] / element.volume_value
            end
            # species after
            for post_element in outneighbors(graph, graph[:name][element.target])
                post_element_index = findall(
                    x ->
                        x.target == props(graph, post_element)[:name] .&&
                        x.source == element.target,
                    elements,
                )[1]
                post_flow_value = elements[post_element_index].flow_value
                du[element_index] -=
                    post_flow_value * u[element_index] / element.volume_value
            end
            # outflow element & outflow -> out
        else #!element.inflow & !startswith("T_", element.source)
            # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
            # or
            # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
            # species before
            for pre_element in inneighbors(graph, graph[:name][element.source])
                pre_element_index = findall(
                    x ->
                        x.source == props(graph, pre_element)[:name] .&&
                        x.target == element.source,
                    elements,
                )[1]
                pre_flow_value = elements[pre_element_index].flow_value
                du[element_index] +=
                    pre_flow_value * u[pre_element_index] / element.volume_value
            end
            # species after
            for _ in outneighbors(graph, graph[:name][element.target])
                du[element_index] -=
                    element.flow_value * u[element_index] / element.volume_value
            end
        end
        # inflow -> terminal node
        if startswith("T_", element.target)
            # dT/dt = QA * A / VT - QV * T / VT
            # species before
            for pre_element in inneighbors(graph, graph[:name][element.target])
                pre_element_index = findall(
                    x ->
                        x.source == props(graph, pre_element)[:name] .&&
                        x.target == element.target,
                    elements,
                )[1]
                pre_flow_value = elements[pre_element_index].flow_value
                du[element_index] +=
                    pre_flow_value * u[pre_element_index] / element.volume_value
            end
            # species after
            for post_element in outneighbors(graph, graph[:name][element.target])
                post_element_index = findall(
                    x ->
                        x.target == props(graph, post_element)[:name] .&&
                        x.source == element.target,
                    elements,
                )[1]
                post_flow_value = elements[post_element_index].flow_value
                du[element_index] -=
                    post_flow_value * u[element_index] / element.volume_value
            end
            # terminal node -> outflow
        elseif startswith("T_", element.source)
            # dV/dt = QV * T / VV - QV * V / VV
            # species before
            for pre_element in inneighbors(graph, graph[:name][element.target])
                pre_element_index = findall(
                    x ->
                        x.source == props(graph, pre_element)[:name] .&&
                        x.target == element.target,
                    elements,
                )[1]
                du[element_index] +=
                    element.flow_value * u[pre_element_index] / element.volume_value
            end
            # species after
            for _ in outneighbors(graph, graph[:name][element.target])
                du[element_index] -=
                    element.flow_value * u[element_index] / element.volume_value
            end
        end
    end
end


# import u0, vessel_tree, p

u0 = [element.initial for element in elements]


tspan = (0.0, 100.0)
prob = ODEProblem(vessel_tree!, u0, tspan)
sol = solve(prob)

header = ["t"]
append!(header, [element.id for element in elements])

df = DataFrame(sol)
rename!(df, Symbol.(header))

CSV.write(
    "/home/mariia/liver_vascularisation/liver_vascular_tree/results/julia_vessel_trees/Rectangle_quad/Rectangle_quad_10/graphs/julia_simulations.csv",
    df,
)

function collect_edges_metadata(; edges_df::DataFrame)::Vector{edge_options}

    # gplot(graph)

    element_suffix::String = ""
    element_id::String = ""

    flow_value::Float64 = 0.0
    volume_value::Float64 = 0.0

    elements = Vector{edge_options}(nothing, length(nrows(edges_df)) + n_terminals)
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
        push!(
            elements,
            edge_options(
                id = element_id,
                source = row.source,
                target = row.target,
                volume_id = "V$(element_id)",
                volume_value = volume_value,
                flow_id = "Q$(element_id)",
                flow_value = flow_value,
                terminal = row.terminal,
                start = row.start,
                inflow = row.is_inflow,
                initial = initial,
            ),
        )
        if startswith("T_", row.target)
            push!(
                elements,
                edge_options(
                    row.target,
                    "",
                    "",
                    volume_terminal,
                    0,
                    true,
                    false,
                    false,
                    initial,
                ),
            )
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

    function f_dxdt(
        dx::Vector{Float64},
        x::Vector{Float64},
        p::Vector{Float64},
        t::Float64,
    )::Nothing
        for (element_index, element) in enumerate(elements)
            # initialize equation
            dx[element_index]::Float64 = 0.0
            # marginal inflow element
            if element.inflow & occursin("marginal", element.id)
                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
                # species after
                for _ in outneighbors(graph, graph[:name][element.target])
                    dx[element_index] -=
                        p[length(elements)+element_index] * x[element_index] /
                        p[element_index]
                end
                # in -> inflow & inflow element
            elseif element.inflow & !occursin("marginal", element.id)
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                # or
                # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                # species before
                for pre_element in inneighbors(graph, graph[:name][element.source])
                    pre_element_index = findall(
                        x ->
                            x.source == props(graph, pre_element)[:name] .&&
                            x.target == element.source,
                        elements,
                    )[1]
                    dx[element_index] +=
                        p[length(elements)+element_index] * x[pre_element_index] /
                        p[element_index]
                end
                # species after
                for post_element in outneighbors(graph, graph[:name][element.target])
                    post_element_index = findall(
                        x ->
                            x.target == props(graph, post_element)[:name] .&&
                            x.source == element.target,
                        elements,
                    )[1]
                    dx[element_index] -=
                        p[length(elements)+post_element_index] * x[element_index] /
                        p[element_index]
                end
                # outflow element & outflow -> out
            else
                !element.inflow & !startswith("T_", element.source)
                # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
                # or
                # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
                # species before
                for pre_element in inneighbors(graph, graph[:name][element.source])
                    pre_element_index = findall(
                        x ->
                            x.source == props(graph, pre_element)[:name] .&&
                            x.target == element.source,
                        elements,
                    )[1]
                    dx[element_index] +=
                        p[length(elements)+pre_element_index] * x[pre_element_index] /
                        p[element_index]
                end
                # species after
                for _ in outneighbors(graph, graph[:name][element.target])
                    dx[element_index] -=
                        p[length(elements)+element_index] * x[element_index] /
                        p[element_index]
                end
            end
            # inflow -> terminal node
            if startswith("T_", element.target)
                # dT/dt = QA * A / VT - QV * T / VT
                # species before
                for pre_element in inneighbors(graph, graph[:name][element.target])
                    pre_element_index = findall(
                        x ->
                            x.source == props(graph, pre_element)[:name] .&&
                            x.target == element.target,
                        elements,
                    )[1]
                    dx[element_index] +=
                        p[length(elements)+pre_element_index] * x[pre_element_index] /
                        p[element_index]
                end
                # species after
                for post_element in outneighbors(graph, graph[:name][element.target])
                    post_element_index = findall(
                        x ->
                            x.target == props(graph, post_element)[:name] .&&
                            x.source == element.target,
                        elements,
                    )[1]
                    dx[element_index] -=
                        p[length(elements)+post_element_index] * x[element_index] /
                        p[element_index]
                end
                # terminal node -> outflow
            elseif startswith("T_", element.source)
                # dV/dt = QV * T / VV - QV * V / VV
                # species before
                for pre_element in inneighbors(graph, graph[:name][element.target])
                    pre_element_index = findall(
                        x ->
                            x.source == props(graph, pre_element)[:name] .&&
                            x.target == element.target,
                        elements,
                    )[1]
                    dx[element_index] +=
                        p[length(elements)+element_index] * x[pre_element_index] /
                        p[element_index]
                end
                # species after
                for _ in outneighbors(graph, graph[:name][element.target])
                    dx[element_index] -=
                        p[length(elements)+element_index] * x[element_index] /
                        p[element_index]
                end
            end
        end
    end
    return f_dxdt!
end
