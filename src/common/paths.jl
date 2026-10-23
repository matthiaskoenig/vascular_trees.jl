module Paths
    export PROJECT_ROOT, JULIA_RESULTS_DIR, BENCHMARKING_RESULTS_PATH

    # repository root (paths.jl lives in src/common/)
    const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
    # absolute, so that the results do not depend on the working directory of the script
    const JULIA_RESULTS_DIR = joinpath(PROJECT_ROOT, "results", "julia_vessel_trees")
    const BENCHMARKING_RESULTS_PATH = joinpath(JULIA_RESULTS_DIR, "jrunning_times.csv")
end
