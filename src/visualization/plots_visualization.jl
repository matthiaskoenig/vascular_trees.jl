module Plots

include("../utils.jl")
import .Utils: BENCHMARKING_RESULTS_PATH

using CSV
using DataFrames
using CairoMakie
using GLMakie
using AlgebraOfGraphics
using DataFramesMeta
using Statistics

function __init__()
    df = CSV.read(BENCHMARKING_RESULTS_PATH, DataFrame)
    plot_run_times_trees(df)
end

function plot_run_times_trees(df)
    df = df[df.call_orders .!= "Total time (all iterations)", :]
    df_error = @by df [:tree_id, :n_terminals, :call_orders, :n_species, :tolerance, :saveat, :version] begin
        :time_min_mean = mean(:times_min)
        :time_min_std = std(:times_min)
        :allocated_gbytes_mean = mean(:allocated_gbytes)
        :allocated_gbytes_std = std(:allocated_gbytes)
    end
    replace!.([df_error.time_min_std, df_error.allocated_gbytes_std], NaN => 0.0)
    sort!(df_error, [:tree_id])

    f1 = Figure()

    #First axis
    plot1 =
        data(df_error) * (
            mapping(
                :n_species => "Number of species nodes, [n]",
                :time_min_mean => "Time, [min]",
                :time_min_std,
                color = :version,
                layout = :saveat => "Tolerance",
            ) * (visual(Scatter) + visual(Errorbars) + smooth())
        )
    subfig1 = draw!(f1[1, 1], plot1)

    #Second axis
    plot2 =
        data(df_error) *
        mapping(
            :n_terminals => "Number of terminal nodes, [n]",
            :allocated_gbytes_mean => "Allocated memory, [GB]",
            :allocated_gbytes_std,
            color = :version,
            layout = :saveat => "Tolerance",
        ) *
        (visual(Scatter) + visual(Errorbars) + smooth())
    draw!(f1[2, 1], plot2)

    # Insert the legend
    legend!(
        f1[0, :],
        subfig1;
        position = :top,
        orientation = :horizontal,
        tellheight = true,
        framevisible = false,
    )

    display(f1)

end

end
