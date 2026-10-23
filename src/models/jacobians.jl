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