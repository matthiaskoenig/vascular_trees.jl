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

include("interventions.jl")
using .Interventions: f_intervention

using ..Definitions: TreeDefinitions, TerminalParameters, VascularTreeParameters

using InteractiveUtils, SparseArrays

# vectors preallocation
const inflow_in = zeros(1)
const inflow_out = zeros(2)
const inflow_out_term = zeros(1)

const outflow_in = zeros(2)
const outflow_out = zeros(1)
const outflow_in_margin = zeros(1)
const outflow_out_margin = zeros(1)

const terminal_in = zeros(2)

export jf_dxdt!, jac!

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

    # TODO: u[end] = f_intervention(t) still writes into u. Replacing it with a separate
    #   dose input means changing how the root's equation reads its predecessor.
    # dose: f_intervention is an infusion rate (amount/min), u holds concentrations, so the
    # marginal element (u[end]) gets the concentration of the incoming blood, rate / Q_root
    # (Q_root: flow of the root edge, the successor of the marginal element)
    root_element = view(post_idx, post_off[end-1]:post_off[end]-1)
    u[end] = f_intervention(t) / sum(view(flows, root_element))
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

function jf_dxdt!(du::Array, u::Array, p::TerminalParameters, t::Float64)
    jf_terminal!(du, u, p, t)
end

function jf_terminal!(du::Array, u::Array, p::TerminalParameters, t::Float64)
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

    `u` and `du` are vectors of length `n_rows * n_terminals`: the `(n_rows, n_terminals)`
        block of `p.flow_values` in column-major order, species `(ke, kt)` at the linear index
        `(kt - 1) * n_rows + ke` (same convention as `jac!`). A vector is needed because the
        analytic-Jacobian callback of Sundials.jl only supports vector states.
    """
    flows = p.flow_values
    volume = p.volumes
    n_rows, n_terminals = size(flows)
    @inbounds for kt in 1:n_terminals
        # linear index of the terminal species (row 1 of column kt)
        kterm = (kt - 1) * n_rows + 1
        inflow = 0.0
        for ke in 2:n_rows
            inflow += flows[ke, kt] * u[kterm + ke - 1]
            du[kterm + ke - 1] = 0.0             # inflow species are set by synchronization
        end
        # dT/dt = Σ Q_inflow * Inflow / VT - Q_outflow * T / VT
        du[kterm] = (inflow - flows[1, kt] * u[kterm]) / volume
    end
end

"""
    jac!(J, p::VascularTreeParameters)

Fill `J` with the Jacobian of `jf_dxdt!` of a vascular tree (`jac_inflow!` or `jac_outflow!`).

Every stored entry of `J` is written, the structural zeros included: the solver does not use
the values of `jac_prototype`, only its sparsity pattern (DiffEqBase replaces the prototype
by `similar(jac_prototype)`, whose values are uninitialized), so an entry that is not written
here would hold random memory.
"""
function jac!(J::SparseArrays.SparseMatrixCSC, p::VascularTreeParameters)
    fill!(nonzeros(J), 0.0)
    if p.is_inflow
        jac_inflow!(J, p)
    else
        jac_outflow!(J, p)
    end
end

"""
    jac_inflow!(J, p)

Jacobian of `jf_inflow!`, `J[i, j] = ∂(du_i)/∂u_j`.

The equations are linear in `u`, so `J` does not depend on `u` or `t`.
It is filled column by column: for every element `ke`, the column `J[:, ke]` holds the
dependence of `ke`'s own equation (diagonal) and of the equations of its successors
(`J[post, ke]`) on `u[ke]`. An element's equation never depends on its successors, so
`J[ke, post] = 0`, and `J` is lower triangular in flow order.

The marginal element (`0`, dose `u[end]`) is a time-dependent input set by
`f_intervention(t)`, not a state: its row and column are zero, its effect enters only
through `∂f/∂t`. Rows of the terminal elements (`3`, `du = 0`) are zero as well.

Only the non-zero entries are assigned; the zeros are written by `jac!`, which calls this
function. Assigning an entry outside the pattern of `J` inserts a new structural entry, which
is how `subsystem_problem` builds the pattern (diagonal + `(post, ke)` pairs) from a diagonal
matrix; on later calls all entries exist and only the values are overwritten.
"""
function jac_inflow!(J::SparseArrays.SparseMatrixCSC, p::VascularTreeParameters)
    flows = p.flow_values
    volumes = p.volume_values
    ODE_groups = p.ODE_groups

    post_idx, post_off = p.successors_idx, p.successors_offset

    @inbounds for (ke, group) in enumerate(ODE_groups)
        post_element = view(post_idx, post_off[ke]:post_off[ke+1]-1)
        # in -> inflow & inflow element not connected to terminal
        if group == 1
            # dA/dt = QA * A_pre / VA - Σ Q_post * A / VA  ->  ∂/∂A = -Σ Q_post / VA
            J[ke, ke] = -sum(view(flows, post_element)) / volumes[ke]
            # successors: dA_post/dt = Q_post * A / V_post - ...  ->  ∂/∂A = Q_post / V_post
            J[post_element, ke] .= view(flows, post_element) ./ view(volumes, post_element)
        # inflow element connected to terminal
        elseif group == 2
            # dA/dt = QA * A_pre / VA - QA * A * n_post / VA  ->  ∂/∂A = -QA * n_post / VA
            # successors are terminals (du = 0), so no off-diagonal entries
            J[ke, ke] = -flows[ke] * length(post_element) / volumes[ke]
        end
    end
end

"""
    jac_outflow!(J, p)

Jacobian of `jf_outflow!`, `J[i, j] = ∂(du_i)/∂u_j`.

The equations are linear in `u`, so `J` does not depend on `u` or `t`.
It is filled column by column: for every element `ke`, the column `J[:, ke]` holds the
dependence of `ke`'s own equation (diagonal) and of the equations of its successors
(`J[post, ke]`) on `u[ke]`. An element's equation never depends on its successors, so
`J[ke, post] = 0`, and `J` is lower triangular in flow order.

Terminal elements (`3`) have zero rows (`du = 0`, their values are set by synchronization
with the terminal part), but non-zero columns: preterminal elements read them.

Only the non-zero entries are assigned; the zeros are written by `jac!`, which calls this
function. Assigning an entry outside the pattern of `J` inserts a new structural entry, which
is how `subsystem_problem` builds the pattern (diagonal + `(post, ke)` pairs) from a diagonal
matrix; on later calls all entries exist and only the values are overwritten.
"""
function jac_outflow!(J::SparseArrays.SparseMatrixCSC, p::VascularTreeParameters)
    flows = p.flow_values
    volumes = p.volume_values
    ODE_groups = p.ODE_groups

    pre_idx, pre_off = p.predecessors_idx, p.predecessors_offset
    post_idx, post_off = p.successors_idx, p.successors_offset

    @inbounds for (ke, group) in enumerate(ODE_groups)
        pre_element = view(pre_idx, pre_off[ke]:pre_off[ke+1]-1)
        post_element = view(post_idx, post_off[ke]:post_off[ke+1]-1)
        # marginal (outflow) element
        if group == 0
            # dV/dt = Σ QV_pre * V_pre / VV - Σ QV_pre * V / VV  ->  ∂/∂V = -Σ QV_pre / VV
            J[ke, ke] = -sum(view(flows, pre_element)) / volumes[ke]
        # outflow -> out & outflow element not connected to terminal,
        # or outflow element connected to terminal
        elseif group == 1 || group == 2
            # dV/dt = ... - QV * V * n_post / VV  ->  ∂/∂V = -QV * n_post / VV
            J[ke, ke] = -flows[ke] * length(post_element) / volumes[ke]
            # successors (groups 0, 1): dV_post/dt = QV * V / V_post - ...  ->  ∂/∂V = QV / V_post
            J[post_element, ke] .= flows[ke] ./ view(volumes, post_element)
        # terminal element (du = 0, set by synchronization)
        elseif group == 3
            # successors (preterminal): dV_post/dt = Q_post * T / V_post - ...  ->  ∂/∂T = Q_post / V_post
            J[post_element, ke] .= view(flows, post_element) ./ view(volumes, post_element)
        end
    end
end

"""
    jac!(J, p::TerminalParameters)

Jacobian of `jf_terminal!`, `J[i, j] = ∂(du_i)/∂u_j`.

`u` is a vector holding the `(n_rows, n_terminals)` block in column-major order, so species
`(ke, kt)` has the linear index `(kt - 1) * n_rows + ke`.
Terminals are independent of each other, so `J` is block diagonal with one
`n_rows * n_rows` block per terminal. In each block only the first row (the terminal
itself) is non-zero: the inflow species (rows 2…) have `du = 0`, they are set by
synchronization with the inflow trees.

Every stored entry of `J` is written, the structural zeros included (see `jac!` for the
trees for the reason). Assigning an entry outside the pattern inserts it, which is how
`subsystem_problem` builds the pattern from a diagonal matrix.
"""
function jac!(J::SparseArrays.SparseMatrixCSC, p::TerminalParameters)
    fill!(nonzeros(J), 0.0)
    flows = p.flow_values
    volume = p.volumes
    n_rows, n_terminals = size(flows)
    @inbounds for kt in 1:n_terminals
        # linear index of the terminal species (row 1 of column kt)
        kterm = (kt - 1) * n_rows + 1
        # dT/dt = Σ Q_inflow * Inflow / VT - Q_outflow * T / VT
        # ∂/∂T = -Q_outflow / VT
        J[kterm, kterm] = -flows[1, kt] / volume
        # ∂/∂Inflow = Q_inflow / VT
        for ke in 2:n_rows
            J[kterm, kterm + ke - 1] = flows[ke, kt] / volume
        end
    end
end


end
