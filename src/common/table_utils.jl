"""
    TableUtils

Small helpers for DataFrames and `.arrow` files that are shared by the trees processing
and the simulation code.
"""
module TableUtils

    export selection_from_df, 
        rows_to_tuples, 
        save_as_arrow, 
        get_extended_vector,
        block_to_df,
        load_arrow

    using DataFrames, 
        Arrow,
        Tables

    """
    selection_from_df(df, (rows, cols)) -> SubDataFrame or SubArray

    View (no copy) of `df[rows, cols]`. `rows` is `:` or a Bool mask, `cols` a column name
    (→ vector view) or a list of names (→ `SubDataFrame`).
    """
    function selection_from_df(df::AbstractDataFrame, (rows, cols)::Tuple)
        return @view df[rows, cols]
    end

    """
    rows_to_tuples(df) -> Vector{Tuple}

    One tuple per row of `df`, e.g. `(source_id, target_id)` or `(x, y, z)`; this is how edges
    and coordinates are stored in the `.arrow` files.
    """
    function rows_to_tuples(df::AbstractDataFrame)
        return Tuple.(Tables.namedtupleiterator(df))
    end

    """
    get_extended_vector(v, n) -> Vector{Union{Missing,T}}

    Copy `v` into a vector of length `n`, padded with `missing` at the end. Used to store
    shorter columns in a table with `n` rows. Throws a `BoundsError` if `length(v) > n`.
    """
    function get_extended_vector(v::AbstractVector{T}, n::Int) where {T}
        out = Vector{Union{Missing,T}}(missing, n)
        copyto!(out, v)
        return out
    end

    """
    block_to_df(block, prefix) -> DataFrame

    Turn one block (`rows * n_terminals`) into a table with one row per terminal and one column
    per block row, named `<prefix>_1`, `<prefix>_2`, …. The block is transposed on its own, so
    the columns keep the element type of the block (`String`, `Float64` or `NTuple{3,Float64}`).
    """
    function block_to_df(block::AbstractMatrix, prefix::String)
        # n_terminals × rows, eltype stays String / Float64 / NTuple{3,Float64}
        return DataFrame(permutedims(block), ["$(prefix)_$k" for k in 1:size(block, 1)])
    end

    """
    save_as_arrow(df, ARROW_PATH)

    Write `df` to `ARROW_PATH` (single-threaded: `ntasks=1`).
    """
    function save_as_arrow(df::DataFrame, ARROW_PATH::String)
        Arrow.write(ARROW_PATH, df, ntasks=1)
    end

    """
    load_arrow(ARROW_PATH[, col_selection]) -> DataFrame

    Load an `.arrow` file as a `DataFrame`, optionally only the columns in `col_selection`.
    """
    function load_arrow(ARROW_PATH::String)::DataFrame
        return DataFrame(Arrow.Table(ARROW_PATH), copycols=false)
    end

    function load_arrow(ARROW_PATH::String, col_selection::Vector{Symbol})::DataFrame
        return DataFrame(Arrow.Table(ARROW_PATH), copycols=false)[
            :,
            col_selection,
        ]
    end
end