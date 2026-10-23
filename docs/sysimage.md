# Sysimage (faster start-up)

- [Build the sysimage](#build-the-sysimage)
- [Use the sysimage](#use-the-sysimage)
  - [Terminal](#terminal)
  - [Visual Studio Code](#visual-studio-code)

See also: [Installation](installation.md) · [Usage](usage.md) · [back to README](../README.md)

Loading `DifferentialEquations`, `Sundials`, `Arrow`, … and compiling the package functions
takes a long time at the start of every Julia session. A custom sysimage built with
[PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl) contains these
packages and the already compiled package functions, so `using VascularTrees` is instant
and the first simulation runs without compilation.

## Build the sysimage

1. Install PackageCompiler once into the global environment (not into the project):
   ```julia
   pkg> activate @v1.13
   pkg> add PackageCompiler
   ```
2. From the repository root, start Julia with the project and build the image:
   ```julia
   julia --project=.
   julia> using PackageCompiler
   julia> create_sysimage(
              [
                  "VascularTrees",
                  # ODE solving
                  "DifferentialEquations", "OrdinaryDiffEq", "Sundials", "SciMLBase", "SparseArrays",
                  # tables / IO
                  "Arrow", "CSV", "DataFrames", "DataFramesMeta", "Tables", "Dictionaries",
                  # utilities used inside the package
                  "Parameters", "Reexport", "ProgressMeter", "TimerOutputs",
              ];
              sysimage_path = "JuliaSysimage.so",
              project = ".",
              precompile_execution_file = "scripts/precompile_project.jl",
          )
   ```
   This takes a while (5-10 min). `JuliaSysimage.so` is written to the repository root and is
   ignored by git.

[scripts/precompile_project.jl](../scripts/precompile_project.jl) is the *precompile workload*:
it processes and simulates the smallest tree (`Rectangle_quad`, 10 nodes). Only the functions
it calls are compiled into the image, so add calls there for every code path you want to
be fast (e.g. other configurations, `benchmark = true`).

**Rebuild the sysimage** after changing code in `src/`, after `Pkg.update`/adding packages,
and after a Julia update. Code baked into the sysimage is *not* updated by Revise; while
developing `src/` intensively, work without the sysimage.

## Use the sysimage

### Terminal

```bash
julia --project -J JuliaSysimage.so
julia --project -J JuliaSysimage.so scripts/run_simulations.jl
```

### Visual Studio Code

Create `.vscode/settings.json` in the repository root (the folder is not tracked by git):

```json
{
    "julia.additionalArgs": ["-J", "${workspaceFolder}/JuliaSysimage.so"]
}
```

- The setting is used only by the REPL started by the Julia extension: `Julia: Start REPL`
  (`Alt+J Alt+O`), or the first `Shift+Enter` in a `.jl` file. A `julia` typed into a normal
  bash terminal ignores it.
- After changing the setting or rebuilding the image, restart the REPL
  (`Julia: Restart REPL`, `Alt+J Alt+R`).
- `${workspaceFolder}` is the first workspace folder, so open the `vascular_trees.jl` folder
  itself in VS Code.
- `julia.useCustomSysimage` is not supported by current versions of the Julia extension.

Check that the sysimage is used:

```julia
julia> unsafe_string(Base.JLOptions().image_file)   # …/vascular_trees.jl/JuliaSysimage.so
```

Note: `Shift+Enter` inside `scripts/run_simulations.jl` or `scripts/build_arrow_graphs.jl`
sends the whole `module … end` block, i.e. starts the full run.
