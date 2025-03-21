module Julia_models

"""
Module with ODE functions that are transport models for the drug(s).
ODEs for inflow trees (arterial and portal) differ from ODEs for
    outflow trees (venous and biliary), so they are written in separate
    functions.
"""

include("../interventions.jl")
using .Interventions: f_intervention

using ..Simulation_Helpers: terminal_inflow, terminal_outflow, terminal_difference

using ...Utils.Definitions: tree_definitions, terminal_parameters, vascular_tree_parameters
using UnsafeArrays

const inflow_in = zeros(1)
const inflow_out = zeros(2)
const inflow_out_term = zeros(1)

const outflow_in = zeros(2)
const outflow_out = zeros(1)
const outflow_in_margin = zeros(1)
const outflow_out_margin = zeros(1)

const terminal_in = zeros(2)

export jf_dxdt!

function jf_dxdt!(du::Vector, u::Vector, p::vascular_tree_parameters, t::Float64)
    if p.is_inflow
        jf_inflow!(du, u, p, t)
    else
        jf_outflow!(du, u, p, t)
    end
end

function jf_inflow!(du, u, p, t)

    flows = p.flow_values
    volumes = p.volume_values
    ODE_groups = p.ODE_groups
    pre_elements = p.pre_elements
    post_elements = p.post_elements
    if p.id == "A"
        u[end] = f_intervention(t)
    end
    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        pre_element = pre_elements[ke]
        post_element = post_elements[ke]
        # in -> inflow & inflow element not connected to terminal
        if group == 1 #(!is_preterminal .&& !is_terminal .&& source_id != 0) 
            # dA/dt = QA * A_pre / VA - Q_post * A / VA;
            # species before
            inflow_in .= flows[ke] .* view(u, pre_element)
            # species after
            inflow_out .= view(flows, post_element) .* u[ke]
            du[ke] = sum(inflow_in) - sum(inflow_out)
            # inflow element connected to terminal
        elseif group == 2 # (is_preterminal) 
            # dA/dt = QA * A_pre / VA - QA * A / VA;
            # species before
            inflow_in .= flows[ke] .* view(u, pre_element)
            du[ke] = sum(inflow_in) - flows[ke] * u[ke] * length(post_element)
            # terminal element
            # elseif group == 3 # (is_terminal)
            #     # THIS IS DIFFERENT THAN IN OTHER MODELS
            #     # dT/dt = QA * A / VT - QA * T / VT
            #     # species before and output
            #     inflow_in .= view(flows, pre_element) .* view(u, pre_element)
            #     # species output
            #     inflow_out_term .= view(flows, pre_element) .* u[ke]
            #     du[ke] = sum(inflow_in) - sum(inflow_out_term)
        end
    end
    du .= du ./ volumes
end

function jf_outflow!(du, u, p, t)

    flows = p.flow_values
    volumes = p.volume_values
    ODE_groups = p.ODE_groups
    pre_elements = p.pre_elements
    post_elements = p.post_elements
    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        pre_element = pre_elements[ke]
        post_element = post_elements[ke]
        # write equations
        # marginal (outflow) element
        if group == 0 #target_id == 0
            # species before
            outflow_in_margin .= view(flows, pre_element) .* view(u, pre_element)
            # output
            outflow_out_margin .= view(flows, pre_element) .* u[ke]
            du[ke] = sum(outflow_in_margin) - sum(outflow_out_margin)
        # outflow -> out & outflow element not connected to terminal
        elseif group == 1 #(!is_preterminal .&& !is_terminal .&& target_id != 0) 
            # dV/dt = QV_pre * V_pre / VV - Q * V / VV;
            # species before
            outflow_in .= view(flows, pre_element) .* view(u, pre_element)
            du[ke] = sum(outflow_in) - flows[ke] * u[ke] * length(post_element)
        # outflow element connected to terminal
        elseif group == 2 #(is_preterminal) 
            # dV/dt = QV * T / VV - QV * V / VV;
            outflow_out .= flows[ke] .* view(u, pre_element)
            # species before and after
            du[ke] =
                sum(outflow_out) -
                flows[ke] * u[ke] * length(post_element)
        elseif group == 3 #(is_terminal) 
            du[ke] = 0
        end
    end
    du .= du ./ volumes
end

function jf_dxdt!(du::Array, u::Array, p::terminal_parameters, t::Float64)
    n_rows = size(u)[1]
    flow_values = uview(p.flow_values)
    du .= 0
    # global terminal_inflow .= view(flow_values, 2:n_rows, :) .* view(u, 2:n_rows, :)
    # global terminal_outflow .= view(flow_values, 1:1, :) .* view(u, 1:1, :)
    # global terminal_difference .= sum(terminal_inflow, dims=1) .- terminal_outflow
    # du[1, :] .= terminal_difference[1, :]
   # @views p.terminal_inflow .= flow_values[2:n_rows, :] .* u[2:n_rows, :]
    p.terminal_inflow .= view(flow_values, 2:n_rows, :) .* view(u, 2:n_rows, :)
    p.terminal_outflow .= view(flow_values, 1:1, :) .* view(u, 1:1, :)
    p.terminal_difference .= sum(p.terminal_inflow, dims=1) .- p.terminal_outflow
    
    du[1, :] .= view(p.terminal_difference, 1, :) ./ p.volumes

end

end
