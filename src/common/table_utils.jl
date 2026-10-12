module TableUtils

    export selection_from_df, 
        create_tuples_from_dfrows, 
        save_as_arrow, 
        get_extended_vector,
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

    function load_arrow(ARROW_PATH::String, col_selection::List)::DataFrame
        # this is a duplicate of a function from get_ode_parameters.jl
        graph_df = DataFrame(Arrow.Table(ARROW_PATH))[
            :,
            col_selection,
        ]
        return graph
    end
end