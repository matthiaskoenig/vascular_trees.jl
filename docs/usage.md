# Usage

- [1. Generate vascular trees](#1-generate-vascular-trees)
- [2. Build arrow graphs](#2-build-arrow-graphs)
- [3. Run simulations](#3-run-simulations)
  - [Options](#options)
  - [Outputs](#outputs)

See also: [Installation](installation.md) · [Project structure](project_structure.md) ·
[Sysimage (faster start-up)](sysimage.md) · [back to README](../README.md)

Before running, make sure that the `results/` folder has the expected layout
(see [Results folder structure](project_structure.md#results-folder-structure)) — the
run scripts do not create any folders themselves.

## 1. Generate vascular trees

Generate the synthetic vascular trees with
[SyntheticVascularTrees.jl](https://gitlab.com/etiennejessen/SyntheticVascularTrees.jl)
(see [Installation](installation.md#install-syntheticvasculartreesjl)). For each vascular
tree (e.g. `A`, `P`, `V`, `B`) three files are expected in the `julia/` subfolder of the
graph directory:

- `<tree>.lg` — edge list (graph skeleton)
- `<tree>_edges.csv` — edge attributes
- `<tree>_nodes.csv` — node attributes

## 2. Build arrow graphs

```bash
julia --project scripts/build_arrow_graphs.jl
```

Processing starts as soon as the file is included. Which trees are processed is set in
the `TreeOptions` struct inside [the script](../scripts/build_arrow_graphs.jl) (`n_nodes`
and `tree_configurations`). The resulting `.arrow` tables are written to the `graphs/`
subfolder of the graph directory.

## 3. Run simulations

```bash
julia --project scripts/run_simulations.jl
```

The simulations start as soon as the file is included, using the options defined in
[the script](../scripts/run_simulations.jl).

### Options

- `TreeOptions` — which tree configurations and numbers of nodes to simulate
  (every combination is run).
- `SimulationOptions` — time span `tspan` [min], number of `steps`, whether to save the
  solutions (`save_simulations`) and whether to benchmark (`benchmark`).
- `SolverOptions` — ODE solver (any
  [SciML algorithm](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/), e.g.
  `Tsit5()` or `CVODE_BDF()`), its tolerances and whether to use the analytic Jacobian
  (`compute_jacobian`, default `true`). The Jacobian is sparse, so with CVODE use a sparse
  linear solver: `CVODE_BDF(linear_solver = :KLU)` (the default solver). KLU requires
  `compute_jacobian = true`; this is checked when `SolverOptions` is constructed.
- `BenchmarkOptions` — number of timed iterations and whether to save the running times
  (used only if `benchmark = true`).
- `flow_scaling_factors` — factors by which the flow is scaled (`1.0` = original flow).

All option structs are defined in [src/common/options.jl](../src/common/options.jl).

### Outputs

In the graph directory:

- `simulations/<subsystem>_simulations_<dt>_dt.arrow` — one file per vascular tree and one
  for the terminal part (if `save_simulations = true` and benchmarking is off). Each table
  has one row per species: a `species_ids` column followed by one column per time point
  (`t_<time>`).
- `results/julia_vessel_trees/jrunning_times.csv` — benchmark running times (if
  `benchmark = true` and `save_running_times = true`).
