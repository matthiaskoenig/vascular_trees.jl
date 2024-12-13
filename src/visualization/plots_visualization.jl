module Plots

    include("../utils.jl")
    import .Utils: BENCHMARKING_RESULTS_PATH

    using CSV
    using DataFrames
    using CairoMakie
    using GLMakie
    using AlgebraOfGraphics

    function __init__()
        plot_running_times(plot_3D=false)
    end

    #allcombinations(v...) = vec(collect(Iterators.product(v...)))

    function plot_running_times(; plot_3D:: Bool)

        df = CSV.read(BENCHMARKING_RESULTS_PATH, DataFrame)
        sort!(df, [:call_orders])

        # find max for time axis, so all plots will have the same time axis grid
        tmax = maximum(df[!, :times_min])
        mmax = maximum(df[!, :times_min])

        if plot_3D
            # ==== 3D plot
            f3 = Figure()
            ax = Axis3(f3[1, 1],
                title = "Running times",
                xlabel = "Execution number, [n]",
                ylabel = "Number of species, [n]",
                zlabel = "Running time, [min]",
                limits = (nothing, nothing, (0, tmax)))

            execution_n = df[!, :n_calls]
            n_species = df[!, :n_species]
            time_min = df[!, :times_min]
            scatter!(ax, execution_n, n_species, time_min; markersize=n_species/100)
            display(f3)
        else
            f2 = Figure()

            #First axis
            plot1 = data(df) * mapping(:call_orders, :times_min, color=:graph_ids) * visual(ScatterLines)
            subfig1 = draw!(f2[1,1], plot1)

            #Second axis
            plot2 = data(df) * mapping(:n_species, :times_min, color=:graph_ids) * visual(ScatterLines)
            subfig2 = draw!(f2[1,2], plot2)

            #Third axis
            plot3 = data(df) * mapping(:call_orders, :allocated_gbytes, color=:graph_ids) * visual(ScatterLines)
            subfig3 = draw!(f2[2,1], plot3)

            #Fourth axis
            plot4 = data(df) * mapping(:n_species, :allocated_gbytes, color=:graph_ids) * visual(ScatterLines)
            subfig4 = draw!(f2[2,2], plot4)

            # Insert the legend
            legend!(
                f2[end+1, 1:2],
                subfig1;
                orientation=:horizontal,
                tellheight=true
            )

            display(f2)

        end

    end
    
end