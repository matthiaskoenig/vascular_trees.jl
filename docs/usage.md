# Usage

- [1. Generate vascular trees](#1-generate-vascular-trees)
- [2. Build arrow graphs](#2-build-arrow-graphs)
- [3. Run simulations](#3-run-simulations)
  - [Options](#options)
  - [Outputs](#outputs)
- [4. Create VTKHDF files (ParaView)](#4-create-vtkhdf-files-paraview)

See also: [Installation](installation.md) · [Project structure](project_structure.md) ·
[Sysimage (faster start-up)](sysimage.md) · [back to README](../README.md)

Only the input files of step 1 must be placed by hand (see
[Results folder structure](project_structure.md#results-folder-structure)); the output
folders of steps 2–4 are created by the run scripts when needed.

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
subfolder of the graph directory (created if needed).

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

- `simulations/<subsystem>[_<f>Q]_simulations_<dt>_dt.arrow` — one file per vascular tree and one
  for the terminal part (if `save_simulations = true` and benchmarking is off). Each table
  has one row per species: a `species_ids` column followed by one column per time point
  (`t_<time>`).
- `results/julia_vessel_trees/jrunning_times.csv` — benchmark running times (if
  `benchmark = true` and `save_running_times = true`).

## 4. Create VTKHDF files (ParaView)

```bash
julia --project scripts/create_vtk.jl
```

Writes the results of step 3 as time-dependent VTKHDF files. Options inside
[the script](../scripts/create_vtk.jl): `TreeOptions`, `flow_scaling_factors` and `dt` — they
select the simulation files, so they must be the values used in step 3
(`dt = (tspan[2] - tspan[1]) / steps`).

Output in the graph directory: `vtkhdf/<tree>[_<f>Q]_<dt>_dt.vtkhdf`, one file per vascular
tree (the terminal part `T` has no geometry of its own). Each file contains the lines (edges)
of the tree with three cell arrays: `Flow` and `Volume` (constant) and `Concentration`
(one value per edge and stored time point).

Open the files in ParaView ≥ 5.12 and color by `Concentration`; the time slider shows the
stored time points. ParaView does not notice when an open file is rewritten — use
File > Reload Files.
