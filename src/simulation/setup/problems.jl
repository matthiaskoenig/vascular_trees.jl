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
        J = jac_prototype(p_sub)
        n_pattern = nnz(J)
        jac!(J, p_sub)
        # jac! must not insert entries outside the pattern (that would be the slow path again)
        @assert nnz(J) == n_pattern "jac_prototype misses entries written by jac!"
        f = ODEFunction(jf_dxdt!; jac = (J, u, p, t) -> jac!(J, p), jac_prototype = J)
    else
        f = jf_dxdt!
    end
    return ODEProblem(f, u0_sub, tspan, p_sub)
end