"""
Functions that set up and run the coupled simulation of the vascular trees and save the results.

Data flow of `run_simulations`:
1. Preparation: ODE parameters, initial values and species ids are loaded for every vascular
   tree (A, P, V, B, ...) and for the terminal part ("T").
2. Integration (`solve_tree!`): time is advanced in steps of `dt`. In every step each tree and
   the terminal part are integrated separately and then coupled through the values at the
   connecting species (see `solve_tree!`).
3. Saving: for each subsystem an ARROW file with one row per species (column `species_ids`) and
   one column per time point (`t_<time>`) is written (only if `sim_options.save_simulations = true`
   and no benchmarking).

Subsystem = vascular tree (e.g. "A") or the terminal part ("T").

The stages are logged with their running times (see `stage`); per-subsystem stages are logged
at `Debug` level, shown with `ENV["JULIA_DEBUG"] = "VascularTrees"`.
"""
module Simulation

    export run_simulations

    using DifferentialEquations, 
        DataFrames, 
        CSV, 
        Dictionaries, 
        TimerOutputs, 
        ProgressMeter, 
        Arrow, 
        SparseArrays,
        Sundials,
        Profile,
        PProf,
        Mmap

    using ..Definitions: FlowDirections, ODE_groups, VascularTreeParameters
    using ..Benchmarking: save_times_as_csv
    using ..TableUtils: save_as_arrow
    using ..LogUtils: stage
    using ..ODEparameters: get_ODE_parameters, get_initial_values
    using ..PharmacokineticModels: jf_dxdt!, jac!, jac_prototype

    # Which trees are inflow / outflow, and the numbers marking the ODE groups of the species
    # (both defined in utils.jl)
    const flow_direction = FlowDirections()
    const groups = ODE_groups()

   include("results/output.jl")
   include("setup/problems.jl")
   include("solver/coupling.jl")
   include("solver/time_loop.jl")
   include("run.jl")

end
