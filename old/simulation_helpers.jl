module Simulation_helpers

export create_simulations, create_benchmarked_simulations

using OrdinaryDiffEq
using DataFrames
using CSV
using TimerOutputs
using ModelingToolkit
using Plots
using Term

include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Definitions: tree_definitions
import .Utils.Benchmarking: save_times_as_csv
import .Utils.Options: model_types

include("../models/julia_from_pygraph.jl")
import .Julia_from_pygraph: get_ODE_components
include("../models/julia_from_jgraph.jl")
import .Julia_from_jgraph: get_ODE_components

# Already specified in utils.jl
m_types::model_types = model_types()
trees::tree_definitions = tree_definitions()

function ODE_solver(;
    ode_system,
    x0,
    tspan,
    tpoints,
    parameter_values,
    sol_options,
    model_type,
)::SciMLBase.ODESolution

    prob = ODEProblem(ode_system, x0, tspan, parameter_values)

    sol = solve(
        prob,
        sol_options.solver, # Rosenbrock23(), # Tsit5(), # CVODE_BDF
        saveat = tpoints,
        dense = false,
        reltol = sol_options.relative_tolerance,
        abstol = sol_options.absolute_tolerance,
    )

    return sol

end

function save_simulations_to_csv(;
    simulations::SciMLBase.ODESolution,
    column_names::Vector{String},
    simulations_path,
)
    # Convert solution to DataFrame
    df = DataFrame(simulations)

    # Write DataFrame to CSV
    header = vcat(["time"], column_names)
    CSV.write(simulations_path, df, header = header)
end

function create_simulations(; t_options, sim_options, sol_options)

    graph_id::String = ""
    file_name::String = ""

    for model_type in t_options.model_types
        if model_type ∈ m_types.templates
            for n_node ∈ t_options.n_nodes, tree_id ∈ t_options.tree_ids
                graph_id = "$(tree_id)_$(n_node)"
                file_name = "$(graph_id)_$(model_type)"
                print("\r...Working with $(file_name)...")
                MODEL_PATH = normpath(
                    joinpath(
                        @__FILE__,
                        "../../..",
                        JULIA_RESULTS_DIR,
                        tree_id,
                        graph_id,
                        "models",
                        "$(file_name).jl",
                    ),
                )
                include(MODEL_PATH) # Load the odes module
                if contains(model_type, "!")
                    f_dxdt = Transport_model.f_dxdt!
                else
                    f_dxdt = Transport_model.f_dxdt
                end
                simulations = ODE_solver(
                    ode_system = f_dxdt,
                    x0 = Transport_model.x0,
                    tspan = sim_options.tspan,
                    tpoints = sim_options.tpoints,
                    parameter_values = Transport_model.p,
                    sol_options = sol_options,
                    model_type = model_type,
                )
                (sim_options.save_simulations) && (save_simulations_to_csv(
                    simulations = simulations,
                    column_names = Transport_model.xids,
                    simulations_path = joinpath(
                        JULIA_RESULTS_DIR,
                        tree_id,
                        graph_id,
                        "simulations",
                        "$(file_name).csv",
                    ),
                ))
            end

        elseif model_type ∈ m_types.julia_model
            MODEL_PATH = normpath(joinpath(@__FILE__, "../../models/Pharmacokinetic_models.jl"))
            include(MODEL_PATH)
            if endswith(model_type, "_python")
                f_dxdt = Pharmacokinetic_models.pyf_dxdt!
                get_ODE_components = Julia_from_pygraph.get_ODE_components
            elseif endswith(model_type, "_julia")
                f_dxdt = Pharmacokinetic_models.jf_dxdt!
                get_ODE_components = Julia_from_jgraph.get_ODE_components
            end
            for n_node ∈ t_options.n_nodes, tree_id ∈ t_options.tree_ids
                graph_id = "$(tree_id)_$(n_node)"
                printstyled(
                    "------------------------------------------------------------------------------------\n";
                    color = 124,
                )
                printstyled("   $(graph_id)_$(model_type)   \n"; color = 9)
                printstyled(
                    "------------------------------------------------------------------------------------\n";
                    color = 124,
                )
                for vessel_tree ∈ trees.vascular_trees[tree_id]
                    x0, p = get_ODE_components(tree_id, n_node, vessel_tree)
                    simulations = ODE_solver(
                        ode_system = f_dxdt,
                        x0 = x0,
                        tspan = sim_options.tspan,
                        tpoints = sim_options.tpoints,
                        parameter_values = p,
                        sol_options = sol_options,
                        model_type = model_type,
                    )
                    display(plot(simulations))
                    if sim_options.save_simulations
                        file_name = "$(graph_id)_$(vessel_tree)_$(model_type)"
                        save_simulations_to_csv(
                            simulations = simulations,
                            column_names = ["$(edge)" for edge in p[3]],
                            simulations_path = joinpath(
                                JULIA_RESULTS_DIR,
                                tree_id,
                                graph_id,
                                "simulations",
                                "$(file_name).csv",
                            ),
                        )
                        #save_simulations_to_csv(simulations=simulations, column_names=["$(elements[1])_$(elements[2])" for elements in eachrow(@view p[:, 1:2])], simulations_path=joinpath(JULIA_RESULTS_DIR, tree_id, graph_id, "simulations", "$(file_name).csv"))
                    end
                end
            end
        end
    end
end

function create_benchmarked_simulations(;
    t_options,
    sim_options,
    bench_options,
    sol_options,
)
    to = TimerOutput()
    graph_id::String = ""
    file_name::String = ""

    pbar = ProgressBar()

    for model_type in t_options.model_types
        if model_type ∈ m_types.templates
            for n_node ∈ t_options.n_nodes, tree_id ∈ t_options.tree_ids
                graph_id = "$(tree_id)_$(n_node)"
                file_name = "$(graph_id)_$(model_type)"
                printstyled(
                    "------------------------------------------------------------------------------------\n";
                    color = 124,
                )
                printstyled("   $(file_name)   \n"; color = 9)
                printstyled(
                    "------------------------------------------------------------------------------------\n";
                    color = 124,
                )

                MODEL_PATH = normpath(
                    joinpath(
                        @__FILE__,
                        "../../..",
                        JULIA_RESULTS_DIR,
                        tree_id,
                        graph_id,
                        "models",
                        "$(file_name).jl",
                    ),
                )
                include(MODEL_PATH) # Load the odes module
                @timeit to "$(graph_id)" begin
                    if contains(model_type, "!")
                        f_dxdt = Transport_model.f_dxdt!
                    else
                        f_dxdt = Transport_model.f_dxdt
                    end
                    @timeit to "$(graph_id)_precompilation" (!contains(model_type, "MT")) &&
                                                            (@invokelatest f_dxdt(
                        zeros(length(Transport_model.x0)),
                        Transport_model.x0,
                        Transport_model.p,
                        0.0,
                    ))
                    for ki in range(1, bench_options.n_iterations, step = 1)
                        @timeit to "$(graph_id)_$ki" ODE_solver(
                            ode_system = f_dxdt,
                            x0 = Transport_model.x0,
                            tspan = sim_options.tspan,
                            tpoints = sim_options.tpoints,
                            parameter_values = Transport_model.p,
                            sol_options = sol_options,
                            model_type = model_type,
                        )
                    end
                end
                n_species = length(Transport_model.x0)
                show(to, sortby = :firstexec)
                (bench_options.save_running_times) && (save_times_as_csv(
                    times = to,
                    n_node = n_node,
                    tree_id = tree_id,
                    n_species = n_species,
                    model_type = model_type,
                    solver_name = sol_options.solver_name,
                    n_term = Int32(0),
                ))
                reset_timer!(to::TimerOutput)
            end

        elseif model_type ∈ m_types.julia_model
            MODEL_PATH = normpath(joinpath(@__FILE__, "../../models/Pharmacokinetic_models.jl"))
            include(MODEL_PATH)
            # if endswith(model_type, "_python")
            #     f_dxdt = Pharmacokinetic_models.pyf_dxdt!
            #     get_ODE_components = Julia_from_pygraph.get_ODE_components
            if endswith(model_type, "_julia")
                f_dxdt = Pharmacokinetic_models.jf_dxdt!
                get_ODE_components = Julia_from_jgraph.get_ODE_components
            end
            for tree_id ∈ t_options.tree_ids, n_node ∈ t_options.n_nodes
                graph_id = "$(tree_id)_$(n_node)"
                printstyled(
                    "------------------------------------------------------------------------------------\n";
                    color = 124,
                )
                printstyled("   $(graph_id)_$(model_type)   \n"; color = 9)
                printstyled(
                    "------------------------------------------------------------------------------------\n";
                    color = 124,
                )
                for vessel_tree ∈ trees.vascular_trees[tree_id]
                    x0, p = get_ODE_components(tree_id, n_node, vessel_tree)
                    printstyled("   ODEs solving started   \n"; color = 9)
                    @timeit to "$(graph_id)_$(vessel_tree)" begin
                        for ki in range(1, bench_options.n_iterations, step = 1)
                            @timeit to "$(graph_id)_$ki" ODE_solver(
                                ode_system = f_dxdt,
                                x0 = x0,
                                tspan = sim_options.tspan,
                                tpoints = sim_options.tpoints,
                                parameter_values = p,
                                sol_options = sol_options,
                                model_type = model_type,
                            )
                        end
                    end
                    (bench_options.save_running_times) && (n_species::Int64 = length(x0))
                    show(to, sortby = :firstexec)
                    (bench_options.save_running_times) && (save_times_as_csv(
                        times = to,
                        n_node = n_node,
                        tree_id = "$(tree_id)_$(vessel_tree)",
                        n_species = n_species,
                        model_type = model_type,
                        solver_name = sol_options.solver_name,
                        n_term = Int32(n_species / 3),
                    ))
                    reset_timer!(to::TimerOutput)
                end
            end
        end
    end
end
end
