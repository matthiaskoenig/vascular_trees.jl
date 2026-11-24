# Presentation: VascularTrees.jl overview

Slides made with [Quarto](https://quarto.org/) (revealjs) and the native Julia engine. The plots
and the example tables are computed from the files in `results/julia_vessel_trees/` each time
the slides are rendered.

## Requirements

- **Quarto ≥ 1.5** (needed for `engine: julia`). Without root rights, install it into `~/.local`:

  ```bash
  VERSION=1.11.9
  curl -LO https://github.com/quarto-dev/quarto-cli/releases/download/v$VERSION/quarto-$VERSION-linux-amd64.tar.gz
  mkdir -p ~/.local/opt ~/.local/bin
  tar xzf quarto-$VERSION-linux-amd64.tar.gz -C ~/.local/opt
  ln -sf ~/.local/opt/quarto-$VERSION/bin/quarto ~/.local/bin/quarto
  quarto --version
  ```

  `~/.local/bin` must be on your `PATH`. On the first render, Quarto installs the Julia package
  that runs the code cells (`QuartoNotebookRunner.jl`) into its own environment.
- **Julia ≥ 1.13**
- **The Julia environment of this folder** (CairoMakie, Arrow, CSV, DataFrames). From the
  repository root, once:

  ```bash
  julia --project=presentation -e 'import Pkg; Pkg.instantiate()'
  ```

- **The results** used by the slides, in `results/julia_vessel_trees/`:
  - `jrunning_times.csv` (benchmark slide)
  - `Rectangle_quad/Rectangle_quad_<n>/graphs/*.arrow` and `.../simulations/*.arrow` for the
    example tables (`n = 10`) and the simulation slides (`n = 100, 10000, 100000, 1000000`)

## Run

From the repository root:

```bash
# live preview: opens the browser and re-renders on every save of a file in this folder
quarto preview presentation/overview.qmd

# or render once to presentation/overview.html
quarto render presentation/overview.qmd
```

A full render takes about 1–2 minutes (Julia start-up and the large simulations). With the
preview running, later renders are faster, because the Julia process is kept alive.

Do not run `quarto render` while `quarto preview` is running for the same file: both use the
same Julia worker and the second one fails with *"the corresponding worker is busy"*. With the
preview running, saving the file is enough.

In VS Code, the **Quarto** extension adds a *Preview* button to `.qmd` files (restart VS Code
after installing Quarto, so that it finds `~/.local/bin/quarto`).

## During the talk

| Key | Action |
|:--|:--|
| `→` / `Space` | next slide or next step (the boxes on the have → how → need slides appear one by one) |
| `F` | full screen |
| `S` | speaker view with the notes and the next slide |
| `O` / `Esc` | overview of all slides |
| `E` | print view (then print to PDF from the browser) |

Slide 2 (*Contents*) links to every section.

The rendered `overview.html` needs the `overview_files/`, `fonts/` and `styles.css` next to it
to be shown elsewhere; copy the whole folder.

## Files

| File | Content |
|:--|:--|
| `overview.qmd` | the slides |
| `plots.jl` | plot and table functions, included in the first (hidden) code cell |
| `styles.css` | layout: have → how → need boxes, cards, contents, matrix diagrams, Roboto font |
| `theme.scss` | revealjs theme variables (fonts of text and headings) |
| `fonts/` | Roboto font files, used by the slides and by the Makie plots (works offline) |
| `Project.toml`, `Manifest.toml` | Julia environment of the slides |

## Changing the plots

The functions are defined in `plots.jl`:

- **Simulation of another graph:** copy one of the *Simulation* slides and change the call

  ```julia
  plot_simulation("Rectangle_quad_PVL_ligated", 1000; dt = 0.02, flow_factor = 1.0,
                  zoom = (6.8, 8.0), max_terminals = 200, size = (1100, 640))
  ```

  `zoom = nothing` removes the zoom panels; `max_terminals` limits the number of terminal curves
  for large trees. `available_simulations()` lists the saved simulations.
- **Benchmark versions:** `plot_benchmark(versions = nothing)` shows the highest `version` in
  `jrunning_times.csv`; pass e.g. `versions = [1, 6, 7]` to compare several.
- **Example tables:** `arrow_example(path; rows, cols)` shows rows and columns of any `.arrow`
  file (path relative to `results/julia_vessel_trees/`).

New section slides (blue): a heading with one `#`, e.g.
`# Simulation results {#results background-color="#2b6cb0"}`. The `{#id}` is what the contents
slide links to (`[text](#results)`).

## Known issues

- Printing to PDF with headless Chrome can leave plots and equations empty; print from the
  normal browser window (`E`, then *Print → Save as PDF*) and check the result.
- The 10⁶ simulation slide shows results of a `T.arrow` built without matching (zero outflow of
  the terminals); it has a warning box until the graph is rebuilt and simulated again.
