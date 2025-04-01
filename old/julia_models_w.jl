module Pharmacokinetic_models

include("../interventions.jl")
using .Interventions: f_intervention

export jf_dxdt!

function jf_dxdt!(dx, x, p, t)

    is_inflow = p[1]
    flows = p[2]
    volumes = p[3]
    ODE_groups = p[4]
    pre_elements = p[5]
    post_elements = p[6]

    if is_inflow == 1
        jf_inflow!(dx, x, flows, volumes, ODE_groups, pre_elements, post_elements, t)
    else
        jf_outflow!(dx, x, flows, volumes, ODE_groups, pre_elements, post_elements, t)
    end
end

function jf_inflow!(dx, x, flows, volumes, ODE_groups, pre_elements, post_elements, t)
    x[length(x)] = f_intervention(t)
    dx .= 0.0
    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        element_flow = flows[ke]
        pre_element = pre_elements[ke]
        post_element = post_elements[ke]
        # in -> inflow & inflow element not connected to terminal
        if group == 1 #(!is_preterminal .&& !is_terminal .&& source_id != 0) 
            # dA/dt = QA * A_pre / VA - Q_post * A / VA;
            # species before and after
            dx[ke] +=
                sum(element_flow .* x[pre_element]) - sum(flows[post_element] .* x[ke])
            # inflow element connected to terminal
        elseif group == 2 # (is_preterminal) 
            # dA/dt = QA * A_pre / VA - QA * A / VA;
            # species before and after
            dx[ke] +=
                sum(element_flow .* x[pre_element]) -
                element_flow * x[ke] * length(post_element)
            # terminal element
        elseif group == 3 # (is_terminal)
            # THIS IS DIFFERENT THAN IN OTHER MODELS
            # dT/dt = QA * A / VT - QA * T / VT
            # species before and output
            dx[ke] +=
                sum(flows[pre_element] .* x[pre_element] .- flows[pre_element] .* x[ke])
        end
    end
    dx .= dx ./ volumes
end

function jf_outflow!(dx, x, flows, volumes, ODE_groups, pre_elements, post_elements, t)
    dx .= 0.0
    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        element_flow = flows[ke]
        pre_element = pre_elements[ke]
        post_element = post_elements[ke]

        # write equations
        # marginal (outflow) element
        if group == 0 #target_id == 0
            # species before
            dx[ke] +=
                sum(flows[pre_element] .* x[pre_element] .- flows[pre_element] .* x[ke])
            # outflow -> out & outflow element not connected to terminal
        elseif group == 1 #(!is_preterminal .&& !is_terminal .&& target_id != 0) 
            # dV/dt = QV_pre * V_pre / VV - Q * V / VV;
            dx[ke] +=
                sum(flows[pre_element] .* x[pre_element]) -
                element_flow .* x[ke] .* length(post_element)
            # outflow element connected to terminal
        elseif group == 2 #(is_preterminal) 
            # dV/dt = QV * T / VV - QV * V / VV;
            dx[ke] +=
                sum(element_flow .* x[pre_element]) -
                element_flow .* x[ke] .* length(post_element)
        end
    end
    dx .= dx ./ volumes
end
end
