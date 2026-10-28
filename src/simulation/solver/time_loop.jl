"""
    solve_tree!(solutions, u0_terminal, p_terminal, u0, p, sim_options, sol_options,
                synch_idxs, additional_sol_options)

Integrate the coupled system of vascular trees and terminal part over `sim_options.tspan`.

Every time step of length `dt`:
1. the values at the start of the step are stored in `solutions` (with adaptive saving only
if they changed enough, see [`integrate!`](@ref));
2. each tree and the terminal part are integrated separately from their current values;
3. the coupling values are exchanged for the next step: the preterminal species of the
inflow trees (A, P) become the inflow of the terminal part, and the first row
of the terminal part (its outflow) becomes the values of the terminal species of
the outflow trees (V, B).

`solutions`, `u0` and `u0_terminal` are modified in place. Returns the stored time points;
`solutions[id][:, k]` is the state at `times[k]`, columns after `length(times)` are unused.

If `sol_options.compute_jacobian`, the trees and the terminal part are solved with the
analytic Jacobian `jac!`.

# Arguments
- `solutions`: for every subsystem a `(n_species, steps+1)` matrix, one column per stored
time point.
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
    
    tmin = sim_options.tspan[1]
    tmax = sim_options.tspan[2]
    dt = sim_options.dt
    # integrate over the whole time span; the loop advances it in steps of dt
    # (tmax + 2dt as margin so that float drift in t never steps past tf)
    tspan = (tmin, tmax + 2dt)

    # no assignments to the outer names inside the closure (they would be boxed)
    integrators, integrator_terminal = stage("setting up integrators") do
        (
            # setup integrators for tree problems
            dictionary(
                id => init_subsystem(id, u0[id], p[id], tspan, sol_options, integrator_options)
                for id in keys(u0)
            ),
            # setup integrator for terminal node problem
            init_subsystem("T", u0_terminal, p_terminal, tspan, sol_options, integrator_options),
        )
    end

    ids = collect(keys(integrators))
    tree_bindings = [(integ = integrators[id], sol = solutions[id], sidx = synch_idxs[id],
                    is_outflow = id in flow_direction.outflow_trees,
                    is_inflow = id in flow_direction.inflow_trees) for id in ids]
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
        enabled = sim_options.show_progress && !sim_options.benchmark
    )
    save_criteria = (abstol = sim_options.save_abstol, reltol = sim_options.save_reltol,
                    max_gap = sim_options.save_max_gap)
    return stage("integration ($(sim_options.steps) steps)") do
        integrate!(solutions["T"], integrator_terminal, size(p_terminal.flow_values, 1),
                inflow_rows, tree_bindings, tmin, dt, save_criteria, progress)
    end

end

"""
    init_subsystem(id, u0_sub, p_sub, tspan, sol_options, integrator_options)

Build the problem of the subsystem `id` (see [`subsystem_problem`](@ref)) and initialize its
integrator. Both steps are logged as `Debug` stages (see [`stage`](@ref)), with the number of
species and, with the analytic Jacobian, the number of stored Jacobian entries.
"""
function init_subsystem(id, u0_sub, p_sub, tspan, sol_options, integrator_options)
    problem = stage("problem $id ($(length(u0_sub)) species)"; debug = true) do
        subsystem_problem(u0_sub, p_sub, tspan, sol_options.compute_jacobian)
    end
    J = problem.f.jac_prototype
    J === nothing || @debug "Jacobian $id: $(nnz(J)) stored entries"
    return stage("init $id"; debug = true) do
        init(problem, sol_options.solver; integrator_options...)
    end
end

"""
    integrate!(solution_terminal, integrator_terminal, n_rows, inflow_rows, tree_bindings,
                tmin, dt, save_criteria, progress) -> times

Time loop of [`solve_tree!`](@ref): store the state of every subsystem at the current time
point (`solution_terminal` and `b.sol` of each tree binding, one column per stored time point),
then advance all subsystems by `dt` and exchange the coupling values. The number of time points
is `size(solution_terminal, 2)`; after the last one is stored, no further step is integrated.

Adaptive saving (`save_criteria = (abstol, reltol, max_gap)`, see `SimulationOptions`): a time
point is stored only if it is the first or the last one, if `max_gap` steps passed since the
last stored one, or if any species of any subsystem changed by more than the tolerance since
the last stored time point (see [`state_changed`](@ref)). All subsystems are stored at the same
time points. With `abstol = reltol = 0` every time point is stored.

Returns the stored time points; the stored states are the first `length(times)` columns.

The marginal element of an inflow tree (`u[end]`, the dose) is not integrated (`du = 0`), and
`jf_inflow!` sets it only in the solver's work vector, so the state would keep its initial
value. It is therefore set to `f_intervention(t)` before every time point is checked and stored.

The terminal state is a vector holding the `(n_rows, n_terminals)` block in column-major
order (see `jf_terminal!`), so row `ki` of the block is the strided range `ki:n_rows:end`.
"""
function integrate!(solution_terminal, integrator_terminal, n_rows, inflow_rows, tree_bindings,
                    tmin, dt, save_criteria, progress)
    n_time_points = size(solution_terminal, 2)
    n_terminal_state = length(integrator_terminal.u)
    save_all = save_criteria.abstol == 0 && save_criteria.reltol == 0
    times = Float64[]
    n_saved = 0
    kl_last_saved = 0
    for kl in 1:n_time_points
        # dose of the inflow trees at the kl-th time point, so that it is stored and seen by
        # the adaptive saving (the solver does not read it back: du = 0, jf_inflow! overwrites it)
        for b in tree_bindings
            b.is_inflow && (b.integ.u[end] = f_intervention(tmin + (kl - 1) * dt))
        end
        # store numerical solutions (state at the kl-th time point) in the next free column
        if save_all || kl == 1 || kl == n_time_points ||
                kl - kl_last_saved >= save_criteria.max_gap ||
                state_changed(tree_bindings, integrator_terminal, solution_terminal, n_saved,
                    save_criteria.abstol, save_criteria.reltol)
            n_saved += 1
            for b in tree_bindings
                copyto!(view(b.sol, :, n_saved), b.integ.u)
            end
            copyto!(view(solution_terminal, :, n_saved), integrator_terminal.u)
            push!(times, tmin + (kl - 1) * dt)
            kl_last_saved = kl
        end
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

        # showvalues as a function: evaluated only when the bar is redrawn
        next!(progress; showvalues = () -> [(:t, tmin + kl * dt), (:stored, n_saved)])
    end
    return times
end

"""
    state_changed(tree_bindings, integrator_terminal, solution_terminal, n_saved, abstol, reltol)

`true` if the current state of any tree or of the terminal part differs from the last stored
column `n_saved` by more than the tolerance, see [`exceeds_tolerance`](@ref).
"""
function state_changed(tree_bindings, integrator_terminal, solution_terminal, n_saved,
                       abstol, reltol)
    for b in tree_bindings
        exceeds_tolerance(b.integ.u, view(b.sol, :, n_saved), abstol, reltol) && return true
    end
    return exceeds_tolerance(integrator_terminal.u, view(solution_terminal, :, n_saved),
                             abstol, reltol)
end

"""
    exceeds_tolerance(u, u_saved, abstol, reltol)

`true` if `|u[i] - u_saved[i]| > abstol + reltol * |u_saved[i]|` for any species `i`
(stops at the first one).
"""
function exceeds_tolerance(u, u_saved, abstol, reltol)
    @inbounds for i in eachindex(u, u_saved)
        abs(u[i] - u_saved[i]) > abstol + reltol * abs(u_saved[i]) && return true
    end
    return false
end

# u was changed from outside the solver. CVODE integrates from its own copy of the state,
# so it is reinitialized from u (Sundials.jl ignores derivative_discontinuity! without callbacks)
reinit_state!(integ::Sundials.CVODEIntegrator) =
    Sundials.CVodeReInit(integ.mem, integ.t, integ.u_nvec)
reinit_state!(integ) = derivative_discontinuity!(integ, true)