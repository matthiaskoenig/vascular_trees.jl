# dev/

Development material that is **not part of the `VascularTrees` package**: it is not loaded by
`using VascularTrees`, not tested, and may not run with the current package version.

- `old/` — legacy code, incl. the benchmarking plots (Makie / AlgebraOfGraphics) in
  `old/visualization/` that are not yet refactored into `src/`.
- `examples/` — early prototype scripts (ODE models of small trees, vectorized and symbolic
  variants, speed-up experiments).
- `sbml/` — experiments with importing SBML models (`.xml`) into Julia.
- `Speedup.md`, `TODO.md` — working notes.
