#=
One-off conversion of the benchmark results (results/julia_vessel_trees/jrunning_times.csv).

The old `n_species` column was the estimate `n_terminals * 9`. It is replaced by the two
counts written by `save_times_as_csv` now (see `count_species` in src/simulation/run.jl):
- `n_states`: states integrated by the solvers;
- `n_distinct_species`: physically distinct species.

Every tree of a Rectangle graph has `2 * n_terminals` edges that are not terminal self-edges
(`n_terminals` preterminal, `n_terminals - 1` inner, 1 marginal), checked on the .arrow files
of all Rectangle_quad graphs. With `n_trees` trees, `n_inflow` of them inflow trees:
- n_distinct_species = (2 * n_trees + 1) * n_terminals
- n_states           = (2 * n_trees + 1 + n_inflow) * n_terminals
(the terminal block has n_inflow + 1 rows, rows 2.. are copies of inflow preterminal species).

The result is written to a new file, the original is not changed.

Run: `julia --project dev/recalculate_species_counts.jl`
=#
using CSV
using DataFrames

const RESULTS_DIR = joinpath(@__DIR__, "..", "results", "julia_vessel_trees")
const INPUT_PATH  = joinpath(RESULTS_DIR, "jrunning_times.csv")
const OUTPUT_PATH = joinpath(RESULTS_DIR, "jrunning_times_recalculated.csv")

# (n_trees, n_inflow) of the configuration in `tree_id`, see TreeDefinitions in definitions.jl;
# Rectangle_quad_A is a run with tree A only
function tree_counts(tree_id::AbstractString)
    startswith(tree_id, "Rectangle_quad_A") && return (n_trees = 1, n_inflow = 1)
    startswith(tree_id, "Rectangle_quad")   && return (n_trees = 4, n_inflow = 2)
    startswith(tree_id, "Rectangle_trio")   && return (n_trees = 3, n_inflow = 2)
    error("Unknown tree configuration: $tree_id")
end

function main()
    # rows appended by the older writer have no tolerance/saveat/step/version (→ missing)
    df = CSV.read(INPUT_PATH, DataFrame)
    counts = tree_counts.(df.tree_id)
    n_trees = getproperty.(counts, :n_trees)
    n_inflow = getproperty.(counts, :n_inflow)

    # the old values must be the distinct count (n_terminals * 9 for quad, * 3 for one tree)
    n_distinct = (2 .* n_trees .+ 1) .* df.n_terminals
    mismatch = findall(df.n_species .!= n_distinct)
    isempty(mismatch) || @warn "Old n_species differs from the distinct count in rows $mismatch"

    # new columns at the position of n_species, the column order of save_times_as_csv
    insertcols!(df, :n_species,
        :n_states => (2 .* n_trees .+ 1 .+ n_inflow) .* df.n_terminals,
        :n_distinct_species => n_distinct)
    select!(df, Not(:n_species))

    CSV.write(OUTPUT_PATH, df)
    @info "Written $(nrow(df)) rows to $OUTPUT_PATH"
end

main()
