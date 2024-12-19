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
        plot_running_times(plot_3D=false)
    end

    function plot_running_times(; plot_3D::Bool)

        df = CSV.read(BENCHMARKING_RESULTS_PATH, DataFrame)
        
        df_error = @by df [:graph_ids, :call_orders, :n_species, :model_type, :algorithm] begin
            :time_min_mean = mean(:times_min)
            :time_min_std = std(:times_min)
            :allocated_gbytes_mean = mean(:allocated_gbytes)
            :allocated_gbytes_std = std(:allocated_gbytes)
        end
        replace!.([df_error.time_min_std, df_error.allocated_gbytes_std], NaN => 0.0)
        sort!(df_error, [:n_species, :model_type])
        

        if !plot_3D

            f2 = Figure()

            #grouped_df = groupby(df_error, :call_orders)
            #First axis
            plot1 = data(df_error) * (
                mapping(:n_species => "Number of species, [n]", 
                :time_min_mean => "Time, [min]", 
                :time_min_std, 
                color=:model_type,
                layout=:call_orders => "Solver name") * 
                (visual(Scatter) + visual(Errorbars) + smooth())
            ) 
            subfig1=draw!(f2[1,1], plot1)

            #Second axis
            plot2 = data(df_error) * mapping(
                :n_species => "Number of species, [n]", 
                :allocated_gbytes_mean => "Allocated memory, [GB]", 
                :allocated_gbytes_std, 
                color=:model_type,
                layout=:call_orders => "Solver name") * 
                (visual(Scatter) + visual(Errorbars) + smooth())
            draw!(f2[2,1], plot2)

            #Label(f2[1, :], "Time with f_dxdt precompilation, first call", fontsize=10, color=:magenta)
            #Label(f2[4, :], "Time without f_dxdt precompilation, first call", fontsize=10, color=:magenta)

            #Label(f2[0, :], "Benchmarking", fontsize=20, color=:magenta)

            # Insert the legend
            legend!(
                f2[0, :],
                subfig1;
                position=:top,
                orientation=:horizontal,
                tellheight=true,
                framevisible=false
            )

            display(f2)
        
        else

            plot3_1 = data(df_error) * (
                mapping(:n_species => "Number of species, [n]", 
                :allocated_gbytes_mean => "Allocated memory, [GB]", 
                :time_min_mean => "Time, [min]", 
                color=:call_orders => "Call order")
            ) 
            plot3 = draw(plot3_1, 
                axis=((type=Axis3, title="Benchmarking", titlesize=20)))

            display(plot3)
        end

    end
    
end