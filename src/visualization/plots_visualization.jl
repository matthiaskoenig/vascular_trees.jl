module Plots

    include("../utils.jl")
    import .Utils: JULIA_RESULTS_DIR

    using CSV
    using Glob
    using DataFrames
    using CairoMakie
    using GLMakie

    function __init__()
        plot_running_times(plot_3D=false)
    end

    function plot_running_times(; plot_3D:: Bool)

        files = glob(joinpath(JULIA_RESULTS_DIR, "running_times_*.csv"))
        dfs = CSV.read.(files, DataFrame)

        # find max for time axis, so all plots will have the same time axis grid
        tmax = 0
        for df ∈ dfs
            tmax_local = maximum(df[!, :times_min])
            if tmax < tmax_local
                tmax=tmax_local
            end
        end

        if plot_3D
            function plot_3D_figure()
                # ==== 3D plot
                f3 = Figure()
                ax = Axis3(f3[1, 1],
                    title = "Running times",
                    xlabel = "Execution number, [n]",
                    ylabel = "Number of species, [n]",
                    zlabel = "Running time, [min]",
                    limits = (nothing, nothing, (0, tmax)))

                for df ∈ dfs
                    execution_n = df[!, :n_calls]
                    n_species = df[!, :n_species]
                    time_min = df[!, :times_min]
                    scatter!(ax, execution_n, n_species, time_min; markersize=n_species/100)
                end
                display(f3)

            end
        else
            function plot_2D()
                # ==== 2D plot
                f2 = Figure()
                kwargs = (; axis = (; limits = (nothing, (0, tmax))))
                kf = 1
                i = 1
                while kf <= length(dfs)
                    j = 1
                    while j < 4 

                        execution_n = dfs[kf][!, :n_calls]
                        n_species = dfs[kf][!, :n_species]
                        time_min = dfs[kf][!, :times_min]
                        scatter(f2[i, j], execution_n, time_min; kwargs..., markersize=n_species/100)
                        j += 1
                        kf += 1
                        if kf > length(dfs)
                            break
                        end

                    end
                    i += 1
                end
                display(f2)
            end
        end

    end
    
end