# TODO

Collected from `src/` and `scripts/`.

## src/trees_processing/terminal_nodes.jl

- [ ] Flow affiliations are wrong. ([L20](src/trees_processing/terminal_nodes.jl#L20))
- [ ] `prepare_tnodes_info` assumes that every tree lists its preterminal edges in the same order, so that the k-th preterminal edge of each tree belongs to terminal k (column k of the output). This is not checked. If it does not hold, values of different terminals get mixed up (possible cause of the wrong flow affiliations). ([L24](src/trees_processing/terminal_nodes.jl#L24))
- [ ] Taken from the "V" tree only; assumes every configuration has a "V" tree. ([L142](src/trees_processing/terminal_nodes.jl#L142))

## scripts/run_simulations.jl

- [ ] Run `@timeit` macros only when needed; currently it works that way, but it is dirty. ([L34](scripts/run_simulations.jl#L34))

