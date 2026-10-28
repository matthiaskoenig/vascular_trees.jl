"""
    TableUtils

Small helpers for DataFrames and `.arrow` files that are shared by the trees processing
and the simulation code: selecting and reshaping tables, loading and saving `.arrow` files, and
`decompose_lists`, which turns a list column (e.g. `pre_elements`) into flat
`(flat_values, offsets)` vectors.
"""
module TableUtils

    export selection_from_df, 
        rows_to_tuples, 
        save_as_arrow, 
        get_extended_vector,
        block_to_df,
        load_arrow,
        decompose_lists

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

    """
    decompose_lists(lists) -> (flat_values, offsets)

    Store a vector of variable-length lists (e.g. `pre_elements` / `post_elements`, one list of
    neighbour indices per edge) as two flat vectors instead of one small vector per list:

    - `flat_values`: all lists concatenated, in order.
    - `offsets`: `length(lists) + 1` start positions; list `k` is
      `flat_values[offsets[k] : offsets[k+1]-1]`. Empty lists give equal neighbouring offsets.

    Example: `[[5], [], [7, 8]]` → `flat_values = [5, 7, 8]`, `offsets = [1, 2, 2, 4]`.

    `lists` can be any vector of vectors, including an Arrow list column, so the per-list vectors
    are never materialized.
    """
    function decompose_lists(lists::AbstractVector{<:AbstractVector{<:Integer}}, first_offset::Int)
        offsets = get_offsets(lists, first_offset)
        flat_values = flatten_lists(lists, offsets[end] - 1)
        return flat_values, offsets
    end
    function decompose_lists(lists::AbstractVector{<:Tuple{<:Integer, <:Integer}}, first_offset::Int)
        vectors_vec = [collect(edge) for edge in lists]
        offsets = get_offsets(vectors_vec, first_offset)
        flat_values = flatten_lists(vectors_vec, offsets[end] - first_offset)
        return flat_values, offsets
    end

    """
    get_offsets(lists) -> Vector{Int}

    Start position of each list in the concatenated vector, plus one entry past the end:
    `offsets[1] = 1`, `offsets[k+1] = offsets[k] + length(lists[k])`. See [`decompose_lists`](@ref).
    """
    function get_offsets(lists::AbstractVector{<:AbstractVector{<:Integer}}, first_offset::Int)::Vector{Int}
        n_lists = length(lists)
        offsets = Vector{Int}(undef, n_lists + 1)
        offsets[1] = first_offset
        @inbounds for k in 1:n_lists
            offsets[k+1] = offsets[k] + length(lists[k])
        end
        return offsets
    end

    """
    flatten_lists(lists, n_values) -> Vector{Int}

    All values of `lists` concatenated in order. `n_values` must be the total number of values
    (`offsets[end] - 1` from [`get_offsets`](@ref)).
    """
    function flatten_lists(lists::AbstractVector{<:AbstractVector{<:Integer}}, n_values::Int)::Vector{Int}
        flat_values = Vector{Int}(undef, n_values)
        position = 1
        @inbounds for list in lists, value in list
            flat_values[position] = value
            position += 1
        end
        return flat_values
    end

end