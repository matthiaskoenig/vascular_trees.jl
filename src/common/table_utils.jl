module TableUtils

    export selection_from_df, 
        create_tuples_from_dfrows, 
        save_as_arrow, 
        get_extended_vector,
        block_to_df,
        load_arrow

    using DataFrames, 
        DataFramesMeta, 
        Arrow

    function selection_from_df(
        df::AbstractDataFrame,
        conditions::Tuple{Union{Colon,BitVector}, Union{Vector{Symbol}, Symbol}},
    )::Union{SubDataFrame, SubArray}
        return @view df[conditions...]
    end

    function create_tuples_from_dfrows(df::AbstractDataFrame)
        return Tuple.(Tables.namedtupleiterator(df))
    end

    function save_as_arrow(df::DataFrame, PATH::String)
        Arrow.write(PATH, df, ntasks=1)
    end

    function get_extended_vector(vector_to_extend, to_which_length_extend::Int)
        return [vector_to_extend; fill(missing, to_which_length_extend - length(vector_to_extend))]
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

    function load_arrow(ARROW_PATH::String, col_selection::Vector{String})::DataFrame
        # this is a duplicate of a function from get_ode_parameters.jl
        graph_df = DataFrame(Arrow.Table(ARROW_PATH))[
            :,
            col_selection,
        ]
        return graph
    end
end