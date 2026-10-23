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
    SimulationOptions(; tspan, steps, save_simulations, benchmark)

    Time settings of the simulation and what to do with the results.

    # Fields
    - `tspan`: `(t_start, t_end)` [min].
    - `steps`: number of time steps; the solution is stored at `steps + 1` time points.
    - `dt`: time step [min], derived from `tspan` and `steps` (DO NOT PASS).
    - `save_simulations`: save the solutions as `.arrow` (only when `benchmark = false`).
    - `benchmark`: run the simulation `n_iterations` times and measure the running time
    instead of saving the solutions, see [`BenchmarkOptions`](@ref).
    """
    @with_kw struct SimulationOptions
        tspan::Tuple{Float64,Float64}
        steps::Int
        dt::Float64 = (tspan[2] - tspan[1]) / steps
        save_simulations::Bool
        benchmark::Bool
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