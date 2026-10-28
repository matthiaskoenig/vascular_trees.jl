# vectors preallocation
const inflow_in = zeros(1)
const inflow_out = zeros(2)
const inflow_out_term = zeros(1)

const outflow_in = zeros(2)
const outflow_out = zeros(1)
const outflow_in_margin = zeros(1)
const outflow_out_margin = zeros(1)

function jf_dxdt!(du::Vector, u::Vector, p::VascularTreeParameters, t::Float64)
    if p.is_inflow
        jf_inflow!(du, u, p, t)
    else
        jf_outflow!(du, u, p, t)
    end
end

function jf_inflow!(du::Vector, u::Vector, p::VascularTreeParameters, t::Float64)
    """
    Function with pharmacokinetic ODE based model for the drug(s) transport along
        liver's vascular inflow (ex. arterial or portal).
    Actors (variables): species in vascular inflow and terminal species.

    Note: Change in terminal species in this function is equal to zero! Here we needed
        them only for right calculation of the preterminal species' change rate. This
        is because right hand side of the equation for the preterminal species includes
        number of outflow ways, where these species should go. So terminal species can not
        be excluded from the inflow's file. So, their rate of change can be any number,
        that will not have any influence. For simplicity, it is equal to zero.
    """
    flows = p.flow_values
    volumes = p.volume_values
    ODE_groups = p.ODE_groups

    pre_idx, pre_off = p.predecessors_idx, p.predecessors_offset
    post_idx, post_off = p.successors_idx, p.successors_offset

    u[end] = f_intervention(t)
    # if p.id == "A"
    #     u[end] = f_intervention(t)
    # end
    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        pre_element = view(pre_idx, pre_off[ke]:pre_off[ke+1]-1)
        post_element = view(post_idx, post_off[ke]:post_off[ke+1]-1)
        # marginal element (dose): set by f_intervention, not integrated
        if group == 0
            du[ke] = 0
        # in -> inflow & inflow element not connected to terminal
        elseif group == 1 #(!is_preterminal .&& !is_terminal .&& source_id != 0)
            # dA/dt = QA * A_pre / VA - Q_post * A / VA;
            # species before
            inflow_in .= flows[ke] .* view(u, pre_element)
            # species after
            inflow_out .= view(flows, post_element) .* u[ke]
            # full equation
            du[ke] = sum(inflow_in) - sum(inflow_out)
        # inflow element connected to terminal
        elseif group == 2 # (is_preterminal) 
            # dA/dt = QA * A_pre / VA - QA * A / VA;
            # species before
            inflow_in .= flows[ke] .* view(u, pre_element)
            # full equation
            du[ke] = sum(inflow_in) - flows[ke] * u[ke] * length(post_element)
        # terminal element
        # for the four vascular trees system is not needed
        elseif group == 3 # (is_terminal)
            du[ke] = 0
        #     # THIS IS DIFFERENT THAN IN OTHER MODELS
        #     # dT/dt = QA * A / VT - QA * T / VT
        #     # species before and output
        #     inflow_in .= view(flows, pre_element) .* view(u, pre_element)
        #     # species output
        #     inflow_out_term .= view(flows, pre_element) .* u[ke]
        #     # full equation
        #     du[ke] = sum(inflow_in) - sum(inflow_out_term)
        end
    end
    du .= du ./ volumes
end

function jf_outflow!(du::Vector, u::Vector, p::VascularTreeParameters, t::Float64)
    """
    Function with pharmacokinetic ODE based model for the drug(s) transport along
        liver's vascular outflow tree (ex. venous or biliary).
    Actors (variables): species in vascular outflow and terminal species.

    Note: Change in terminal species in this function is equal to zero! However,
        in comparison to the jf_inflow! function, here we need them not for the
        calculation of the number of inflow or outflow ways. Here the right hand
        side of the equation for the preterminal species includes amount of the
        terminal species (includes elements from u vector that represent terminal
        species). Nevertheless, outflow file does not contain all elements for
        calculation of the rate change of terminal species, while this is a separate
        function and terminal and outflow parts are being synchronized using terminal
        elements from u vector. So, their rate of change must be a constant.
    """
    flows = p.flow_values
    volumes = p.volume_values
    ODE_groups = p.ODE_groups

    pre_idx, pre_off = p.predecessors_idx, p.predecessors_offset
    post_idx, post_off = p.successors_idx, p.successors_offset

    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        pre_element = view(pre_idx, pre_off[ke]:pre_off[ke+1]-1)
        post_element = view(post_idx, post_off[ke]:post_off[ke+1]-1)
        # write equations
        # marginal (outflow) element
        if group == 0 #target_id == 0
            # dV/dt = QVpre * V_pre / VV - Q * V / VV;
            # species before
            outflow_in_margin .= view(flows, pre_element) .* view(u, pre_element)
            # output
            outflow_out_margin .= view(flows, pre_element) .* u[ke]
            # full equation
            du[ke] = sum(outflow_in_margin) - sum(outflow_out_margin)
        # outflow -> out & outflow element not connected to terminal
        elseif group == 1 #(!is_preterminal .&& !is_terminal .&& target_id != 0) 
            # dV/dt = QV_pre * V_pre / VV - Q * V / VV;
            # species before
            outflow_in .= view(flows, pre_element) .* view(u, pre_element)
            # full equation
            du[ke] = sum(outflow_in) - flows[ke] * u[ke] * length(post_element)
        # outflow element connected to terminal
        elseif group == 2 #(is_preterminal) 
            # dV/dt = QV * T / VV - QV * V / VV;
            outflow_out .= flows[ke] .* view(u, pre_element)
            # full equation
            du[ke] =
                sum(outflow_out) -
                flows[ke] * u[ke] * length(post_element)
        # terminal element
        elseif group == 3 #(is_terminal) 
            du[ke] = 0
        end
    end
    du .= du ./ volumes
end