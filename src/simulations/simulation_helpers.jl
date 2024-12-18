module Simulation_helpers

    export create_simulations, create_benchmarked_simulations

    using OrdinaryDiffEq
    using DataFrames
    using CSV
    using TimerOutputs

    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Benchmarking: save_times_as_csv
    import .Utils.Options: solver_options

    function ODE_solver(; ode_system, x0, tspan, tpoints, parameter_values, sol_options)::SciMLBase.ODESolution
        
        prob = ODEProblem(
            ode_system, 
            x0,
            tspan,
            parameter_values)

        sol = solve(
            prob, 
            sol_options.solver, # Rosenbrock23(), # Tsit5(), # CVODE_BDF
            saveat=tpoints,
            dense=false,
            reltol=sol_options.relative_tolerance, 
            abstol=sol_options.absolute_tolerance,
            dt=sol_options.dt)

        return sol

    end

    function save_simulations_to_csv(; simulations::SciMLBase.ODESolution, column_names::Vector{String}, simulations_path)
        # Step 3: Convert solution to DataFrame
        df = DataFrame(simulations)

        # Step 4: Write DataFrame to CSV
        header = vcat(["time"], column_names)
        CSV.write(simulations_path, df, header=header)
    end

    function create_simulations(; g_options, sim_options)
        
        graph_id::String = ""
        file_name::String = ""
        for n_node ∈ g_options.n_nodes, tree_id ∈ g_options.tree_ids
                graph_id = "$(tree_id)_$(n_node)"
                file_name = "$(graph_id)_$(g_options.model_type)"
                MODEL_PATH = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "models", "$(file_name).jl"))
                include(MODEL_PATH) # Load the odes module
                if contains(g_options.model_type, "!")
                    simulations = ODE_solver(ode_system=Transport_model.f_dxdt!, x0=Transport_model.x0, tspan=sim_options.tspan, tpoints=sim_options.tpoints, parameter_values=Transport_model.p, sol_options=sol_options)
                else
                    simulations = ODE_solver(ode_system=Transport_model.f_dxdt, x0=Transport_model.x0, tspan=sim_options.tspan, tpoints=sim_options.tpoints, parameter_values=Transport_model.p, sol_options=sol_options)
                end

                (sim_options.save_simulations) && (save_simulations_to_csv(simulations=simulations, column_names=Transport_model.xids, simulations_path=joinpath(JULIA_RESULTS_DIR, tree_id, graph_id, "simulations", "$(file_name).csv")))
        end
    end

    function create_benchmarked_simulations(; g_options, sim_options, bench_options, sol_options)
        to = TimerOutput()
        graph_id::String = ""
        file_name::String = ""
        
        for model_type in g_options.model_types, n_node ∈ g_options.n_nodes, tree_id ∈ g_options.tree_ids
                n_species::Array{Int32} = []
                graph_id = "$(tree_id)_$(n_node)"
                file_name = "$(graph_id)_$(model_type)"
                print("\r...Working with $(file_name)...")
                
                MODEL_PATH = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "models", "$(file_name).jl"))
                include(MODEL_PATH) # Load the odes module\
                @timeit to "$(graph_id)" begin
                    if contains(model_type, "!")
                        @invokelatest Transport_model.f_dxdt!(zeros(size(Transport_model.x0)...), Transport_model.x0, Transport_model.p, 0.0)
                        for ki in range(1, bench_options.n_iterations, step=1)
                                @timeit to "$(graph_id)_$ki" ODE_solver(ode_system=Transport_model.f_dxdt!, x0=Transport_model.x0, tspan=sim_options.tspan, tpoints=sim_options.tpoints, parameter_values=Transport_model.p, sol_options=sol_options)
                        end
                    else
                        @invokelatest Transport_model.f_dxdt(zeros(size(Transport_model.x0)...), Transport_model.x0, Transport_model.p, 0.0)
                        for ki in range(1, bench_options.n_iterations, step=1)
                                @timeit to "$(graph_id)_$ki" ODE_solver(ode_system=Transport_model.f_dxdt, x0=Transport_model.x0, tspan=sim_options.tspan, tpoints=sim_options.tpoints, parameter_values=Transport_model.p, sol_options=sol_options)
                        end
                    end
                end
                (bench_options.save_running_times) && (push!(n_species, length(Transport_model.xids)))
                show(to, sortby=:firstexec)
                #(sol_options.krylov) && (model_type="$(model_type)_krylov")
                (bench_options.save_running_times) && (save_times_as_csv(times=to, n_species=n_species, model_type=model_type, solver_name=sol_options.solver_name))
                reset_timer!(to::TimerOutput)
        end
    end
end