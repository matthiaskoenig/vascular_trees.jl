module Simulation_helpers

    export create_simulations, create_benchmarked_simulations

    using OrdinaryDiffEq
    using DataFrames
    using CSV
    using TimerOutputs

    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR
    import .Utils.Benchmarking: save_times_as_csv

    absolute_tolerance = 1e-6
    relative_tolerance = 1e-6

    function ODE_solver(; ode_system, x0, tspan, tpoints, parameter_values)::SciMLBase.ODESolution

        prob = ODEProblem(
            ode_system, 
            x0,
            tspan, 
            parameter_values)

        sol = solve(
            prob, 
            Tsit5(), # Rosenbrock23(), # Tsit5(), # CVODE_BDF
            saveat=tpoints,
            dense=false,
            reltol=relative_tolerance, 
            abstol=absolute_tolerance)

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
        for n_node ∈ g_options.n_nodes
            for tree_id ∈ g_options.tree_ids

                graph_id = "$(tree_id)_$(n_node)"
                file_name = "$(graph_id)_$(g_options.file_suffix)"
                MODEL_PATH = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "models", "$(file_name).jl"))
                include(MODEL_PATH) # Load the odes module
                simulations = ODE_solver(ode_system=Transport_model.f_dxdt, x0=Transport_model.x0, tspan=sim_options.tspan, tpoints=sim_options.tpoints, parameter_values=Transport_model.p)

                if sim_options.save_simulations
                    save_simulations_to_csv(simulations=simulations, column_names=Transport_model.xids, simulations_path=joinpath(JULIA_RESULTS_DIR, tree_id, graph_id, "simulations", "$(file_name).csv"))
                end
            end
        end
    end

    function create_benchmarked_simulations(; g_options, sim_options, bench_options)
        to = TimerOutput()
        kp::Int16 = 1
        graph_id::String = ""
        file_name::String = ""
        n_species::Array{Int32} = []
        for n_node ∈ g_options.n_nodes
            for tree_id ∈ g_options.tree_ids
                for _ in range(1, bench_options.n_iterations, step=1)

                    graph_id = "$(tree_id)_$(n_node)"
                    file_name = "$(graph_id)_$(g_options.file_suffix)"
                    MODEL_PATH = normpath(joinpath(@__FILE__, "../../.." , JULIA_RESULTS_DIR, tree_id, graph_id, "models", "$(file_name).jl"))
                    @timeit to "$(graph_id)_$(kp)" begin
                        include(MODEL_PATH) # Load the odes module
                        simulations = ODE_solver(ode_system=Transport_model.f_dxdt, x0=Transport_model.x0, tspan=sim_options.tspan, tpoints=sim_options.tpoints, parameter_values=Transport_model.p)

                        if sim_options.save_simulations
                            save_simulations_to_csv(simulations=simulations, column_names=Transport_model.xids, simulations_path=joinpath(JULIA_RESULTS_DIR, tree_id, graph_id, "simulations", "$(file_name).csv"))
                        end
                        if bench_options.save_running_times
                            push!(n_species, length(Transport_model.xids))
                        end
                    end

                    kp = kp + 1
                end
            end
        end
        show(to)

        if bench_options.save_running_times
            save_times_as_csv(times=to, n_species=n_species)
        end
    end

end