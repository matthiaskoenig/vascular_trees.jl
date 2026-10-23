const terminal_in = zeros(2)

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