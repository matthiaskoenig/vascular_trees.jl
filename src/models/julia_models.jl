module Julia_models

"""
Module with ODE functions that are transport models for the drug(s).
ODEs for inflow trees (arterial and portal) differ from ODEs for
    outflow trees (venous and biliary), so they are written in separate
    functions.
"""

include("../interventions.jl")
using .Interventions: f_intervention
using InteractiveUtils

const inflow_in = zeros(1)
const inflow_out = zeros(2)
const inflow_out_term = zeros(1)

const outflow_in = zeros(2)
const outflow_in_margin = zeros(1)
const outflow_out_margin = zeros(1)

const terminal_in = zeros(1)

export jf_dxdt!

function jf_dxdt!(dx, x, p, t)
    is_inflow = p[1]
    if is_inflow
        jf_inflow!(dx, x, p, t)
    else
        jf_outflow!(dx, x, p, t)
    end
end

function jf_inflow!(dx, x, p, t)
    
    flows = p[2]
    volumes = p[3]
    ODE_groups = p[4]
    pre_elements = p[5]
    post_elements = p[6]
    x[end] = f_intervention(t)
    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        pre_element = pre_elements[ke]
        post_element = post_elements[ke]
        # in -> inflow & inflow element not connected to terminal
        if group == 1 #(!is_preterminal .&& !is_terminal .&& source_id != 0) 
            # dA/dt = QA * A_pre / VA - Q_post * A / VA;
            # species before
            inflow_in .= flows[ke] .* view(x, pre_element)
            # species after
            inflow_out .= view(flows, post_element) .* x[ke]
            dx[ke] = sum(inflow_in) - sum(inflow_out)
        # inflow element connected to terminal
        elseif group == 2 # (is_preterminal) 
            # dA/dt = QA * A_pre / VA - QA * A / VA;
            # species before
            inflow_in .= flows[ke] .* view(x, pre_element)
            dx[ke] = sum(inflow_in) - flows[ke] * x[ke] * length(post_element)
        # terminal element
        elseif group == 3 # (is_terminal)
            # THIS IS DIFFERENT THAN IN OTHER MODELS
            # dT/dt = QA * A / VT - QA * T / VT
            # species before and output
            inflow_in .= view(flows, pre_element) .* view(x, pre_element)
            # species output
            inflow_out_term .= view(flows, pre_element) .* x[ke]
            dx[ke] = sum(inflow_in) - sum(inflow_out_term)
        end
    end
    dx .= dx ./ volumes
end

function jf_outflow!(dx, x, p, t)
    is_inflow = p[1]
    flows = p[2]
    volumes = p[3]
    ODE_groups = p[4]
    pre_elements = p[5]
    post_elements = p[6]
    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        pre_element = pre_elements[ke]
        post_element = post_elements[ke]
        # write equations
        # marginal (outflow) element
        if group == 0 #target_id == 0
            # species before
            outflow_in_margin .= view(flows, pre_element) .* view(x, pre_element)
            # output
            outflow_out_margin .= view(flows, pre_element) .* x[ke]
            dx[ke] = sum(outflow_in_margin) - sum(outflow_out_margin)
        # outflow -> out & outflow element not connected to terminal
        elseif group == 1 #(!is_preterminal .&& !is_terminal .&& target_id != 0) 
            # dV/dt = QV_pre * V_pre / VV - Q * V / VV;
            # species before
            outflow_in .= view(flows, pre_element) .* view(x, pre_element)
            dx[ke] = sum(outflow_in) - flows[ke] * x[ke] * length(post_element)
        # outflow element connected to terminal
        elseif group == 2 #(is_preterminal) 
            # dV/dt = QV * T / VV - QV * V / VV;
            # species before and after
            dx[ke] = sum(flows[ke] .* view(x, pre_element)) - flows[ke] * x[ke] * length(post_element)
        elseif group == 3 #(is_terminal) 
            dx[ke] = f_intervention(t) 
        end
    end
    dx .= dx ./ volumes
end

function jf_terminal!(dx, x, flows, ODE_groups, pre_elements, post_elements, t)
    @inbounds for (ke, group) in enumerate(ODE_groups)
        if group == 3 # terminal
            pre_element = pre_elements[ke]
            post_element = post_elements[ke]
            terminal_in .= view(flows, pre_element) .* view(x, pre_element)
            terminal_out .= view(flows, post_element) .* x[ke]
            dx[ke] = sum(terminal_in) - sum(terminal_out)
        end
    end
end

end
