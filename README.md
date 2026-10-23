# Vascular trees in Julia

`VascularTrees.jl` is a Julia package for processing synthetic vascular networks of the
liver and for building and simulating pharmacokinetic computational models on these
graphs. The vascular trees are generated with
[SyntheticVascularTrees.jl](https://gitlab.com/etiennejessen/SyntheticVascularTrees.jl);
this package converts them into an efficient tabular format and runs coupled ODE
simulations of substance transport along the trees.

The corresponding Python code is located at
https://github.com/matthiaskoenig/vascular_trees.

- **Tree processing:** converts generated trees into `.arrow` tables.
- **Simulation:** coupled pharmacokinetic ODE simulations from the inflow trees (arterial,
  portal) through the liver sinusoids to the outflow trees (venous, biliary), with optional
  benchmarking.
- **Visualization:** benchmarking plots and mapping of results onto the VTK geometry
  (legacy code in `old/`).

More details: [Project structure](docs/project_structure.md).

## Table of contents

- [Installation](#installation)
- [Usage](#usage)
- [Examples / Demos](#examples--demos)
- [Documentation](#documentation)
- [License](#license)
- [Contributors](#contributors)
- [Funding](#funding)

## Installation

Requirements: [Julia](https://julialang.org/) ≥ 1.13 and
[SyntheticVascularTrees.jl](https://gitlab.com/etiennejessen/SyntheticVascularTrees.jl)
for generating the trees.

```bash
git clone https://github.com/matthiaskoenig/vascular_trees.jl.git
cd vascular_trees.jl
julia --project -e "import Pkg; Pkg.instantiate()"
```

Details (SyntheticVascularTrees.jl, environment activation, VS Code, multithreading):
[docs/installation.md](docs/installation.md). For a faster start-up build a sysimage:
[docs/sysimage.md](docs/sysimage.md).

## Usage

```bash
# 1. convert the generated trees into .arrow tables
julia --project scripts/build_arrow_graphs.jl
# 2. run the coupled-tree simulations
julia --project scripts/run_simulations.jl
```

The options (tree configurations, numbers of nodes, solver, benchmarking, ...) are set
inside the scripts. The scripts read from and write into an existing
[`results/` folder layout](docs/project_structure.md#results-folder-structure).
All options and outputs: [docs/usage.md](docs/usage.md).

## Examples / Demos

Minimal example on the smallest tree (`Rectangle_quad`, 10 nodes), the same as
[scripts/precompile_project.jl](scripts/precompile_project.jl). The input files of the tree
must exist in `results/julia_vessel_trees/Rectangle_quad/Rectangle_quad_10/julia/`.

```julia
using VascularTrees
using VascularTrees.ProcessIndividualTrees: process_julia_graph
using VascularTrees.ProcessTerminalNodes: process_terminal_nodes
using VascularTrees.Simulation: run_simulations

trees = TreeDefinitions()
cfg = "Rectangle_quad"
tree_info = TreeStructure(; tree_configuration = cfg, n_node = 10,
                          tree_components = trees.vascular_trees[cfg])

# build the .arrow tables
process_julia_graph(tree_info)
process_terminal_nodes(tree_info)

# simulate 16 min with 800 time steps and save the solutions
run_simulations(
    tree_info,
    SimulationOptions(tspan = (0.0, 16.0), steps = 800,
                      save_simulations = true, benchmark = false),
    SolverOptions(),
    (dense = false, save_everystep = false, progress = false),
    1.0,                 # flow scaling factor
    BenchmarkOptions(),
)
```

## Documentation

- [Installation](docs/installation.md): requirements, SyntheticVascularTrees.jl,
  dependencies, environment activation, multithreading
- [Usage](docs/usage.md): workflow, all options of the run scripts, outputs
- [Project structure](docs/project_structure.md): features, supported tree configurations,
  repository and results folder structure
- [Sysimage](docs/sysimage.md): building and using a sysimage (terminal and VS Code)

## License

- Source Code: [LGPL-3.0](http://opensource.org/licenses/LGPL-3.0)
- Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

This program is free software: you can redistribute it and/or modify it under the terms
of the GNU Lesser General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

## Contributors

- Mariia Myshkina
- Matthias König

## Funding

Mariia Myshkina and Matthias König are supported by the Federal Ministry of Education
and Research (BMBF, Germany) within the ATLAS project (grant number 031L0304B).

Matthias König is supported by the Federal Ministry of Education and Research (BMBF,
Germany) within the research network Systems Medicine of the Liver (**LiSyM**, grant
number 031L0054) and by the German Research Foundation (DFG) within the Research Unit
Programme FOR 5151 "[QuaLiPerF](https://qualiperf.de) (Quantifying Liver
Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach)" by
grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311,
Subproject SimLivA).

© 2024–2026 Mariia Myshkina & Matthias König
