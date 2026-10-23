# Installation

- [Requirements](#requirements)
- [Clone the repository](#clone-the-repository)
- [Install SyntheticVascularTrees.jl](#install-syntheticvasculartreesjl)
- [Install dependencies](#install-dependencies)
- [Environment activation](#environment-activation)
  - [Activate the project](#activate-the-project)
  - [Multithreading](#multithreading)

See also: [Usage](usage.md) · [Sysimage (faster start-up)](sysimage.md) · [back to README](../README.md)

## Requirements

- [Julia](https://julialang.org/) ≥ 1.13
- Optional: [Visual Studio Code](https://code.visualstudio.com/) with the
  [Julia extension](https://www.julia-vscode.org/)

## Clone the repository

```bash
git clone https://github.com/matthiaskoenig/vascular_trees.jl.git
cd vascular_trees.jl
```

## Install SyntheticVascularTrees.jl

This package depends on the
[SyntheticVascularTrees.jl](https://gitlab.com/etiennejessen/SyntheticVascularTrees.jl)
package developed by Etienne Jessen and
[Dominik Schillinger](https://www.bauing.tu-darmstadt.de/fachbereich_bau_umwelt/ueber_den_fachbereich/professoren____innen_am_fb/professorenliste_details_115392.de.jsp)
at the TU Darmstadt for generating the vascular trees.

- clone the repository and open a terminal in `SyntheticVascularTrees.jl`
- start `julia`, enter the package mode with `]`, then:
  `activate .` → `build` → `test`
- choose the Julia environment in Visual Studio Code

## Install dependencies

All dependencies of this package are specified in `Project.toml` / `Manifest.toml` and
are installed with:

```bash
julia --project -e "import Pkg; Pkg.instantiate()"
```

## Environment activation

### Activate the project

Either start Julia with the project flag from the repository root:

```bash
julia --project=.
```

or activate the environment from the Julia package mode:

```julia
julia> ]
pkg> activate .
```

In Visual Studio Code, select the environment via the command palette
(`Julia: Change Current Environment`) or the environment indicator in the status bar.

### Multithreading

To allow Julia to use multiple threads, set the `JULIA_NUM_THREADS` environment
variable, e.g.:

```bash
JULIA_NUM_THREADS=auto julia --project=.
```

or equivalently start Julia with `julia --threads auto`. In Visual Studio Code this can
be configured in `settings.json`:

```json
"julia.NumThreads": "auto"
```
