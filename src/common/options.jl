"""
    Options

Option structs set by the user in the run scripts (run_simulations.jl,
build_arrow_graphs.jl, map_vtk.jl): which trees to use, how to simulate,
which ODE solver to use and whether to benchmark.
"""
module Options
    export TreeOptions, 
        SimulationOptions, 
        BenchmarkOptions, 
        SolverOptions

    using Parameters
    using Sundials
    using OrdinaryDiffEq

    """
    TreeOptions(; n_nodes, tree_configurations)

    Which trees to process / simulate. Every combination of the two lists is run.

    # Fields
    - `n_nodes`: numbers of nodes of the trees, e.g. `[10, 100, 1000]`.
    - `tree_configurations`: configuration names, e.g. `["Rectangle_quad"]`; must be keys of
    `TreeDefinitions().vascular_trees`.
    """
    @with_kw struct TreeOptions
        n_nodes::Vector{Int64}
        tree_configurations::Vector{String}
    end

    """
    SimulationOptions(; tspan, steps, save_simulations, benchmark,
                       save_abstol = 0.0, save_reltol = 0.0, save_max_gap = typemax(Int),
                       show_progress = true)

    Time settings of the simulation and what to do with the results.

    # Fields
    - `tspan`: `(t_start, t_end)` [min].
    - `steps`: number of time steps (coupling steps of length `dt`); the solution can be stored
    at up to `steps + 1` time points.
    - `dt`: time step [min], derived from `tspan` and `steps` (DO NOT PASS).
    - `save_simulations`: save the solutions as `.arrow` (only when `benchmark = false`).
    - `benchmark`: run the simulation `n_iterations` times and measure the running time
    instead of saving the solutions, see [`BenchmarkOptions`](@ref).
    - `save_abstol`, `save_reltol`: adaptive saving. A time point is stored only if a species of
    any subsystem differs from the last stored time point by more than
    `save_abstol + save_reltol * |u_saved|`. If both are `0.0`, every time point is stored.
    - `save_max_gap`: with adaptive saving, a time point is stored at least every
    `save_max_gap` steps. The first and the last time points are always stored.
    - `show_progress`: show a progress bar of the time loop (simulated time, number of stored
    time points); never shown when `benchmark = true`.
    """
    @with_kw struct SimulationOptions
        tspan::Tuple{Float64,Float64}
        steps::Int
        dt::Float64 = (tspan[2] - tspan[1]) / steps
        save_simulations::Bool
        benchmark::Bool
        save_abstol::Float64 = 0.0
        save_reltol::Float64 = 0.0
        save_max_gap::Int = typemax(Int)
        show_progress::Bool = true
        @assert save_abstol >= 0 && save_reltol >= 0 "save_abstol and save_reltol must be >= 0"
        @assert save_max_gap >= 1 "save_max_gap must be >= 1"
    end

    """
        BenchmarkOptions(; save_running_times = false, n_iterations = 1)

    Used only if `SimulationOptions.benchmark = true`.

    # Fields
    - `save_running_times`: write the measured running times to a CSV file.
    - `n_iterations`: how many times every simulation is run and timed.
    """
    @with_kw struct BenchmarkOptions
        save_running_times::Bool = false
        n_iterations::Int = 1
    end

    """
    SolverOptions(; solver = Tsit5(), absolute_tolerance = 1e-8,
                   relative_tolerance = 1e-8, solver_name = string(nameof(typeof(solver))),
                   compute_jacobian = true)

    ODE solver passed to `init` for every tree and the terminal part, see
    https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/.

    # Fields
    - `solver`: any SciML ODE algorithm, e.g. `Tsit5()` or `CVODE_BDF()` (Sundials).
    - `absolute_tolerance`, `relative_tolerance`: `abstol` / `reltol` of the integrator.
    - `solver_name`: name of the solver written to the benchmark results.
    - `compute_jacobian`: build the analytic Jacobian of every tree (see `jac!`). It is sparse,
      so CVODE needs a sparse linear solver: `CVODE_BDF(linear_solver = :KLU)`.
      KLU needs the sparsity pattern of the Jacobian, so a KLU solver requires
      `compute_jacobian = true` (checked on construction).
    """
    # https://docs.sciml.ai/DiffEqDocs/stable/solvers/split_ode_solve/
    @with_kw struct SolverOptions{A<:SciMLBase.AbstractODEAlgorithm}
        solver::A = CVODE_BDF(linear_solver = :KLU) # "Tsit5", "CVODE_BDF", …
        absolute_tolerance::Float64 = 1e-8
        relative_tolerance::Float64 = 1e-8
        solver_name::String = string(nameof(typeof(solver)))
        compute_jacobian::Bool = true
        @assert !uses_klu(solver) || compute_jacobian "linear_solver = :KLU needs compute_jacobian = true"
    end

    # whether the solver is a Sundials algorithm with the sparse KLU linear solver
    uses_klu(::Sundials.SundialsODEAlgorithm{M, L}) where {M, L} = L == :KLU
    uses_klu(::SciMLBase.AbstractODEAlgorithm) = false

end