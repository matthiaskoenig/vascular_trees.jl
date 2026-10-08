module Options
    export tree_options, simulations_options, benchmark_options, solver_options, edge_options
    using Parameters

    using OrdinaryDiffEq
    using Revise # this package must not be in final version
    using Sundials

    @with_kw struct tree_options
        n_nodes::Vector{Int64}
        tree_configurations::Vector{String}
    end

    Base.@kwdef struct simulations_options
        tspan::Tuple{Float64,Float64}
        steps::Integer
        dt::Float64 = tspan[2] / steps
        save_simulations::Bool
        benchmark::Bool
    end

    @with_kw struct benchmark_options
        save_running_times::Bool = false
        n_iterations::Int16 = 1
    end

    # https://docs.sciml.ai/DiffEqDocs/stable/solvers/split_ode_solve/
    @with_kw struct solver_options
        solver = Tsit5()
        absolute_tolerance = 1e-8
        relative_tolerance = 1e-8
        solver_name = "Tsit5" # "Tsit5"
    end

end