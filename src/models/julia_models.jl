module Pharmacokinetic_models
"""
Module with pharmacokinetic ODE based models for the drug(s): transport of the
    substance(s) along the liver's vessel tree (arterial/portal vascular trees ->
    liver sinusoids (terminal part) -> venous/biliary trees).

ODEs for inflow trees (arterial and portal) and for terminal part differ from ODEs for
    outflow trees (venous and biliary), so they are written in separate functions.
Also, variables (initial values and differentials) and parts of parameters structure 
    differ between ODE functions for vascular trees (vectors) and for terminal part 
    (arrays).
"""

include("../interventions.jl")
using .Interventions: f_intervention

# using ..Simulation_Helpers: terminal_inflow, terminal_outflow, terminal_difference

using ...Utils.Definitions: tree_definitions, terminal_parameters, vascular_tree_parameters
using InteractiveUtils

# vectors preallocation
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

function jf_inflow!(du::Vector, u::Vector, p::vascular_tree_parameters, t::Float64)
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
    pre_elements = p.pre_elements
    post_elements = p.post_elements
    u[end] = f_intervention(t)
    # if p.id == "A"
    #     u[end] = f_intervention(t)
    # end
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

function jf_outflow!(du::Vector, u::Vector, p::vascular_tree_parameters, t::Float64)
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
    pre_elements = p.pre_elements
    post_elements = p.post_elements
    @inbounds for (ke, group) in enumerate(ODE_groups)
        # retrieve information for element
        pre_element = pre_elements[ke]
        post_element = post_elements[ke]
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

function jf_dxdt!(du::Array, u::Array, p::terminal_parameters, t::Float64)
    jf_terminal!(du, u, p, t)
end

function jf_terminal!(du::Array, u::Array, p::terminal_parameters, t::Float64)
    """
    Function with pharmacokinetic ODE based model for the drug(s) transport to 
        sinusoid from liver's inflows and from sinusoid.
    Actors (variables): species in vascular inflow and terminal species.

    Note: Change in inflow species in this function is equal to zero! Here the 
        right hand side of the equation for the preterminal species includes amount 
        of the inflow species (includes elements from u array that represent inflow
        species). Nevertheless, terminal part file does not contain all elements for
        calculation of the rate change of preterminal inflow species, while this is 
        a separate function and terminal and inflow parts are being synchronized using 
        preterminal inflow elements from u array. So, their rate of change must be a 
        constant.
    """
    n_ODE_elements = size(u)[1]
    du .= 0
    # # terminal_difference = p.terminal_difference
    # flow_values = p.flow_values
    # sign = 1
    # ke = 1
    # kt = 1
    # @inbounds for flow_value in flow_values
    #     if ke == 1
    #         sign = -1
    #     else
    #         sign = 1
    #     end
    #     du[1, kt] = du[1, kt] + sign * flow_value * u[ke, kt]
    #     ke += 1
    #     if ke == 4
    #         ke = 1
    #         kt += 1
    #     end
    # end
    # du[1, :] ./= p.volumes

    # dT/dt = Q_inflow * Inflow / VT - Q_outflow * T / VT
    # species before
    p.terminal_inflow .= view(p.flow_values, 2:n_ODE_elements, :) .* view(u, 2:n_ODE_elements, :)
    # output
    p.terminal_outflow .= view(p.flow_values, 1:1, :) .* view(u, 1:1, :)
    # full equation
    p.terminal_difference .= sum(p.terminal_inflow, dims=1) .- p.terminal_outflow
    
    du[1, :] = p.terminal_difference ./ p.volumes
end
end
