include("../utils.jl")
import .Utils: JULIA_RESULTS_DIR
import .Utils.Options: graph_options, simulations_options

include("simulation_helpers.jl")
import .Simulation_helpers: ODE_solver

# ============ Specify options
g_options::graph_options = graph_options(
    n_nodes = [10],
    tree_ids = [
        #"Rectangle_quad",
        "Rectangle_trio",
    ],
    model_type = "symbolic", #"vectorized" "symbolic"
)

tspan = (0.0, 10.0 / 60)
sim_options::simulations_options = simulations_options(
    tspan = tspan,
    tpoints = range(tspan[1], stop = tspan[2], length = 1001),
    save_simulations = false,
    benchmark = false,
)

graph_id::String = ""
file_name::String = ""
for n_node ∈ g_options.n_nodes
    for tree_id ∈ g_options.tree_ids

        graph_id = "$(tree_id)_$(n_node)"
        file_name = "$(graph_id)_$(g_options.model_type)"
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
        import .Transport_model: f_dxdt, x0, p
        simulations = ODE_solver(
            ode_system = f_dxdt,
            x0 = x0,
            tspan = sim_options.tspan,
            tpoints = sim_options.tpoints,
            parameter_values = p,
        )

        # if sim_options.save_simulations
        #     save_simulations_to_csv(simulations=simulations, column_names=Transport_model.xids, simulations_path=joinpath(JULIA_RESULTS_DIR, tree_id, graph_id, "simulations", "$(file_name).csv"))
        # end
    end
end
