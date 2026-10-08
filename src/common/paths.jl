module Paths
    export RESULTS_DIR, JULIA_RESULTS_DIR, BENCHMARKING_RESULTS_PATH, MODEL_PATH
    
    RESULTS_DIR::String = "results"
    JULIA_RESULTS_DIR::String = RESULTS_DIR * "/julia_vessel_trees"
    BENCHMARKING_RESULTS_PATH::String = joinpath(JULIA_RESULTS_DIR, "jrunning_times.csv")
    MODEL_PATH::String = "src/models/pharmacokinetic_models.jl"
end