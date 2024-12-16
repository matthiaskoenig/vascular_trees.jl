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
        plot_running_times(plot_3D=true)
    end

    function plot_running_times(; plot_3D::Bool)

        df = CSV.read(BENCHMARKING_RESULTS_PATH, DataFrame)
        
        df_error = @by df [:graph_ids, :call_orders, :n_species, :tree_ids] begin
            :time_min_mean = mean(:times_min)
            :time_min_std = std(:times_min)
            :allocated_gbytes_mean = mean(:allocated_gbytes)
            :allocated_gbytes_std = std(:allocated_gbytes)
        end
        sort!(df_error, [:n_species])
        replace!.([df_error.time_min_std, df_error.allocated_gbytes_std], NaN => 0.0)

        if !plot_3D

            f2 = Figure()

            #First axis
            plot1 = data(df_error) * (
                mapping(:n_species => "Number of species, [n]", 
                :time_min_mean => "Time, [min]", 
                :time_min_std, 
                color=:call_orders => "Call order") * 
                (visual(Scatter) + visual(Errorbars) + smooth())
            ) 
            subfig1 = draw!(f2[3,1], plot1)

            #Second axis
            plot2 = data(df_error) * mapping(
                :n_species => "Number of species, [n]", 
                :allocated_gbytes_mean => "Allocated memory, [GB]", 
                :allocated_gbytes_std, 
                color=:call_orders => "Call order") * 
                (visual(Scatter) + visual(Errorbars) + smooth())
            draw!(f2[3,2], plot2)

            Label(f2[0, :], "Benchmarking", fontsize=20, color=:magenta)

            # Insert the legend
            legend!(
                f2[1, :],
                subfig1;
                position=:top,
                orientation=:horizontal,
                tellheight=true,
                framevisible=false
            )
            
            Label(f2[2, 1], "Running time", fontsize=10)
            Label(f2[2, 2], "Memory allocation", fontsize=10)

            display(f2)
        
        else

            plot1 = data(df_error) * (
                mapping(:n_species => "Number of species, [n]", 
                :allocated_gbytes_mean => "Allocated memory, [GB]", 
                :time_min_mean => "Time, [min]", 
                color=:call_orders => "Call order")
            ) 
            plot3 = draw(plot1, 
                axis=((type=Axis3, title="Benchmarking", titlesize=20)))

            display(plot3)
        end

    end
    
end