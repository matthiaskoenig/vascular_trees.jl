"""
        solve_tree!(solutions, u0_terminal, p_terminal, u0, p, sim_options, sol_options,
                    synch_idxs, additional_sol_options)

    Integrate the coupled system of vascular trees and terminal part over `sim_options.tspan`.

    Every time step of length `dt`:
    1. each tree and the terminal part are integrated separately from their current values;
    2. the values at the start of the step are stored in `solutions`;
    3. the coupling values are exchanged for the next step: the preterminal species of the
    inflow trees (A, P) become the inflow of the terminal part, and the first row
    of the terminal part (its outflow) becomes the values of the terminal species of
    the outflow trees (V, B).

    `solutions`, `u0` and `u0_terminal` are modified in place.

    If `sol_options.compute_jacobian`, the trees and the terminal part are solved with the
    analytic Jacobian `jac!`.

    # Arguments
    - `solutions`: for every subsystem `steps+1` vectors, one per time point, each holding the
    values of all species of that subsystem.
    - `u0_terminal`, `p_terminal`: initial values (a vector, see `jf_terminal!`) and parameters
    of the terminal part.
    - `u0`, `p`: initial values and parameters of the trees, with tree ids as keys.
    - `sim_options`, `sol_options`, `additional_sol_options`: see [`run_simulations`](@ref).
    - `synch_idxs`: indices of the species that connect each tree with the terminal part
    (see [`get_synchronization_indices`](@ref)).
    """
    function solve_tree!(
        solutions,
        u0_terminal,
        p_terminal,
        u0,
        p,
        sim_options,
        sol_options,
        synch_idxs,
        additional_sol_options,
    )

        # collecting integrator arguments together
        integrator_options = (
            reltol = sol_options.relative_tolerance,
            abstol = sol_options.absolute_tolerance,
            save_start = false,
            save_end = false,
            maxiters = (sim_options.steps + 1) * 100_000,
            additional_sol_options...
        )
        
        @info "Running integration"
        tmin = sim_options.tspan[1]
        tmax = sim_options.tspan[2]
        dt = sim_options.dt
        # integrate over the whole time span; the loop advances it in steps of dt
        # (tmax + 2dt as margin so that float drift in t never steps past tf)
        tspan = (tmin, tmax + 2dt)

        # setup integrators for tree problems
        integrators = dictionary(
            id => init(subsystem_problem(u0[id], p[id], tspan, sol_options.compute_jacobian),
                    sol_options.solver; integrator_options...)
            for id in keys(u0)
        )
        # setup integrator for terminal node problem
        problem_terminal =
            subsystem_problem(u0_terminal, p_terminal, tspan, sol_options.compute_jacobian)
        integrator_terminal = init(problem_terminal, sol_options.solver; integrator_options...)

        ids = collect(keys(integrators))
        tree_bindings = [(integ = integrators[id], sol = solutions[id], sidx = synch_idxs[id],
                        is_outflow = id in flow_direction.outflow_trees) for id in ids]
        # terminal rows that receive inflow from an inflow tree: (row index, tree id)
        # computed once, because the assignment of terminal rows to trees does not change
        inflow_rows = [(ki, findfirst(==(tid), ids)) for (ki, tid) in
               enumerate(first.(view(p_terminal.species_ids, :, 1), 1)) if
               tid in flow_direction.inflow_trees]   # binding index instead of String

        progress = Progress(
            Int(sim_options.steps);
            dt = 1.0,
            barglyphs = BarGlyphs('|', '█', ['▁', '▂', '▃', '▄', '▅', '▆', '▇'], ' ', '|'),
            color = :magenta,
            enabled = false
        )
        integrate!(solutions["T"], integrator_terminal, size(p_terminal.flow_values, 1),
                   inflow_rows, tree_bindings, dt, progress)

    end

    """
        integrate!(solution_terminal, integrator_terminal, n_rows, inflow_rows, tree_bindings,
                   dt, progress)

    Time loop of [`solve_tree!`](@ref): store the state of every subsystem at each time point
    (`solution_terminal` and `b.sol` of each tree binding, one vector per time point), then
    advance all subsystems by `dt` and exchange the coupling values. The number of time points
    is `length(solution_terminal)`; after the last one is stored, no further step is integrated.

    The terminal state is a vector holding the `(n_rows, n_terminals)` block in column-major
    order (see `jf_terminal!`), so row `ki` of the block is the strided range `ki:n_rows:end`.
    """
    function integrate!(solution_terminal, integrator_terminal, n_rows, inflow_rows, tree_bindings,
                        dt, progress)
        n_time_points = length(solution_terminal)
        n_terminal_state = length(integrator_terminal.u)
        for kl in 1:n_time_points
            # store numerical solutions (state at the kl-th time point)
            for b in tree_bindings
                copyto!(b.sol[kl], b.integ.u)
            end
            copyto!(solution_terminal[kl], integrator_terminal.u)
            # last time point stored: nothing left to integrate
            kl == n_time_points && break

            # solve current step: advance every subsystem exactly to t + dt
            for b in tree_bindings
                step!(b.integ, dt, true)
            end
            step!(integrator_terminal, dt, true)

            # updating values at the connecting species
            # update in terminal part (inflow from the inflow trees): row ki of the block
            for (ki, bi) in inflow_rows
                b = tree_bindings[bi]
                copyto!(view(integrator_terminal.u, ki:n_rows:n_terminal_state), view(b.integ.u, b.sidx))
            end
            # update for outflow trees (outflow from the terminal part): row 1 of the block
            # reset its cached derivative
            for b in tree_bindings
                if b.is_outflow
                    b.integ.u[b.sidx] .= view(integrator_terminal.u, 1:n_rows:n_terminal_state)
                    reinit_state!(b.integ)
                end
            end
            # reset its cached derivative
            reinit_state!(integrator_terminal)

            next!(progress)
        end
    end

    # u was changed from outside the solver. CVODE integrates from its own copy of the state,
    # so it is reinitialized from u (Sundials.jl ignores derivative_discontinuity! without callbacks)
    reinit_state!(integ::Sundials.CVODEIntegrator) =
        Sundials.CVodeReInit(integ.mem, integ.t, integ.u_nvec)
    reinit_state!(integ) = derivative_discontinuity!(integ, true)

    """
        subsystem_problem(u0_sub, p_sub, tspan, compute_jacobian) -> ODEProblem

    `ODEProblem` of one subsystem (a vascular tree or the terminal part) for
    [`jf_dxdt!`](@ref), used by [`solve_tree!`](@ref). The method of `jf_dxdt!` and `jac!` is
    chosen by the type of `p_sub` (`VascularTreeParameters` or `TerminalParameters`).

    If `compute_jacobian` is `true`, the problem gets the analytic Jacobian [`jac!`](@ref).
    `jac_prototype` only conveys the sparsity pattern: DiffEqBase replaces it by
    `similar(jac_prototype)` (values uninitialized) and the solver fills that copy via `jac`
    before every factorization, so `jac!` writes every stored entry. Here `jac!` is called once
    on a diagonal matrix to build the pattern (its assignments insert the off-diagonal entries).
    Otherwise the solver approximates the Jacobian itself.

    With CVODE, use a sparse linear solver (`CVODE_BDF(linear_solver = :KLU)`), otherwise the
    sparse Jacobian is converted to a dense one.

    # Arguments
    - `u0_sub`, `p_sub`: initial values and parameters of the subsystem. `u0_sub` must be a
      vector (the Jacobian callback of Sundials.jl does not support matrix states); for the
      terminal part this is the column-major block, see `jf_terminal!`.
    - `tspan`: integration time span.
    - `compute_jacobian`: use the analytic Jacobian, `sol_options.compute_jacobian`.
    """
    function subsystem_problem(u0_sub, p_sub, tspan, compute_jacobian::Bool)
        if compute_jacobian
            n = length(u0_sub)
            # the diagonal is stored even where it is zero: CVODE builds I - γJ in place of J,
            # and with KLU the sparsity pattern of J must already contain it
            J = sparse(1:n, 1:n, zeros(n), n, n)
            jac!(J, p_sub)      # builds the pattern
            f = ODEFunction(jf_dxdt!; jac = (J, u, p, t) -> jac!(J, p), jac_prototype = J)
        else
            f = jf_dxdt!
        end
        return ODEProblem(f, u0_sub, tspan, p_sub)
    end