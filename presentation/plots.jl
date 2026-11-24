# Plots for overview.qmd (benchmark, simulations). Included once in the setup chunk.
using Arrow, CSV, DataFrames, CairoMakie, Statistics, Printf, Markdown

const RESULTS_DIR = normpath(joinpath(@__DIR__, "..", "results", "julia_vessel_trees"))

# ================== Common style
# Labels placed in the figure grid must not "tell" their size along the axes they span
# (tellwidth = false for horizontal, tellheight = false for rotated labels): otherwise the
# column / row is shrunk to the size of the label text and the axes become tiny.
const TREE_COLORS = Dict(
    "A" => "#D81818",   # arterial - red
    "P" => "#8E44AD",   # portal - purple
    "V" => "#0C5BB0",   # venous - blue
    "B" => "#6CC04A",   # biliary - light green
    "T" => "#D57326",   # terminal part - orange
)
const SERIES_COLORS = ["#7209b7", "#f72585", "green", "#0C5BB0", "#D57326", "#4cc9f0"]
const LABEL_STYLE = (font = :bold, color = "#27476e", fontsize = 18)
const SUBTITLE_STYLE = (titlefont = :regular, titlecolor = "#004e89", titlesize = 16)
const ZOOM_COLOR = ("#f72585", 0.2)

const ROBOTO = (
    regular = joinpath(@__DIR__, "fonts", "Roboto-Regular.ttf"),
    bold = joinpath(@__DIR__, "fonts", "Roboto-Bold.ttf"),
)
set_theme!(Theme(fonts = ROBOTO, fontsize = 14,
                 Axis = (xgridcolor = "#EEEEEE", ygridcolor = "#EEEEEE",
                         xminorgridvisible = true, yminorgridvisible = true,
                         xminorgridcolor = "#F5F5F5", yminorgridcolor = "#F5F5F5")))

# ================== Benchmark
"""
    plot_benchmark(; versions = nothing, file = "jrunning_times.csv")

Running time and allocated memory per call against the number of species, like
`plot_jrunning_times` (Python): rows = time / memory, columns = solver tolerances, one line
(mean ± std of the calls) per version and solver. `versions = nothing` plots the highest
version in the file; pass e.g. `[1, 6, 7]` to compare several.
"""
function plot_benchmark(; versions = nothing, file = "jrunning_times.csv")
    df = CSV.read(joinpath(RESULTS_DIR, file), DataFrame)
    # columns added later by save_times_as_csv do not exist in files with only old rows
    hasproperty(df, :jacobian) || (df.jacobian = fill(missing, nrow(df)))
    versions = something(versions, [maximum(skipmissing(df.version))])

    calls = filter(r -> r.call_orders == "Call" && coalesce(r.version in versions, false) &&
                        !ismissing(r.tolerance) && !ismissing(r.solver), df)
    # old rows have the Jacobian in the solver name, new rows in the `jacobian` column
    calls.solver_label = [r.solver * (coalesce(r.jacobian, false) ? " + Jacobian" : "")
                          for r in eachrow(calls)]
    stats = combine(groupby(calls, [:version, :solver_label, :tolerance, :n_distinct_species]),
                    [:times_min, :allocated_gbytes] .=> mean,
                    [:times_min, :allocated_gbytes] .=> (x -> length(x) > 1 ? std(x) : 0.0) .=>
                        [:times_min_std, :allocated_gbytes_std])
    sort!(stats, :n_distinct_species)

    tolerances = sort(unique(stats.tolerance); rev = true)
    series = unique(stats[:, [:version, :solver_label]])
    several_solvers(v) = count(==(v), series.version) > 1
    series_label(v, s) = "Version $v" * (several_solvers(v) ? " · $s" : "")

    fig = Figure(size = (1000, 640))
    Label(fig[0, 1:length(tolerances)], "Model Performance"; font = :bold, fontsize = 18, tellwidth = false)
    quantities = [(:times_min, "Time, [min]"), (:allocated_gbytes, "Allocated memory, [GB]")]
    axes = [Axis(fig[row, col]; xscale = log10, yscale = log10,
                 title = row == 1 ? "Tolerance = $(tolerances[col])" : "", SUBTITLE_STYLE...)
            for row in 1:2, col in eachindex(tolerances)]
    for (row, (_, ylabel)) in enumerate(quantities)
        Label(fig[row, 0], ylabel; rotation = π / 2, tellheight = false, LABEL_STYLE...)
        linkyaxes!(axes[row, :]...)
    end
    Label(fig[3, 1:length(tolerances)], "Species, [n]"; tellwidth = false, LABEL_STYLE...)

    xmax = maximum(stats.n_distinct_species) * 2
    for ax in axes
        lines!(ax, [1, 1e9], [1, 1e9]; color = :black, linestyle = :dot, label = "Identity line")
    end
    for (ks, s) in enumerate(eachrow(series)), (col, tol) in enumerate(tolerances)
        d = stats[(stats.version .== s.version) .& (stats.solver_label .== s.solver_label) .&
                  (stats.tolerance .== tol), :]
        isempty(d) && continue
        color = SERIES_COLORS[mod1(ks, length(SERIES_COLORS))]
        for (row, (q, _)) in enumerate(quantities)
            μ, σ = d[!, Symbol(q, "_mean")], d[!, Symbol(q, "_std")]
            ax = axes[row, col]
            # lower error bar must stay > 0 on a log axis
            errorbars!(ax, d.n_distinct_species, μ, min.(σ, 0.9 .* μ), σ; color, whiskerwidth = 6)
            scatterlines!(ax, d.n_distinct_species, μ; color, markersize = 10,
                          label = series_label(s.version, s.solver_label))
        end
    end
    # limits from the data, so that the identity line does not stretch them
    for (row, (q, _)) in enumerate(quantities)
        μ, σ = stats[!, Symbol(q, "_mean")], stats[!, Symbol(q, "_std")]
        lo, hi = minimum(μ .- min.(σ, 0.9 .* μ)), maximum(μ .+ σ)
        for ax in axes[row, :]
            limits!(ax, 10, xmax, lo / 2, hi * 2)
        end
    end
    n_entries = nrow(series) + 1
    Legend(fig[4, 1:length(tolerances)], axes[1, 1]; orientation = :horizontal,
           nbanks = n_entries > 4 ? 2 : 1, tellwidth = false, framecolor = "#004e89",
           unique = true)
    return fig
end

# ================== Example tables
"""
    arrow_example(path; rows, cols) -> Markdown.MD

Markdown table with the rows `rows` (a `⋮` row is put where rows are skipped) and the columns
`cols` of the `.arrow` file `path` (relative to `RESULTS_DIR`), for the slides on the file layout.
"""
function arrow_example(path; rows, cols)
    df = DataFrame(Arrow.Table(joinpath(RESULTS_DIR, path)))
    fmt(x) = x isa Missing ? "missing" :
             x isa AbstractFloat ? @sprintf("%.3g", x) :
             x isa AbstractVector ? "[" * join(x, ", ") * "]" :
             x isa String ? "`$x`" : string(x)
    lines = ["| " * join(("`$c`" for c in cols), " | ") * " |",
             "|" * repeat(":--:|", length(cols))]
    rows = rows isa Function ? rows(nrow(df)) : rows
    for (k, r) in enumerate(rows)
        k > 1 && r != rows[k-1] + 1 && push!(lines, "| " * join(fill("⋮", length(cols)), " | ") * " |")
        push!(lines, "| " * join((fmt(df[r, c]) for c in cols), " | ") * " |")
    end
    return Markdown.parse(join(lines, "\n"))
end

# ================== Simulations
"""
    available_simulations() -> DataFrame

Graphs (configuration, number of nodes) and time steps for which simulations are saved.
"""
function available_simulations()
    found = DataFrame(tree_configuration = String[], n_node = String[], file = String[])
    for cfg in readdir(RESULTS_DIR), graph in (isdir(joinpath(RESULTS_DIR, cfg)) ? readdir(joinpath(RESULTS_DIR, cfg)) : String[])
        sim_dir = joinpath(RESULTS_DIR, cfg, graph, "simulations")
        isdir(sim_dir) || continue
        for f in readdir(sim_dir)
            startswith(f, "T") && endswith(f, ".arrow") &&
                push!(found, (cfg, chopprefix(graph, cfg * "_"), f))
        end
    end
    return found
end

"""
    read_species(path, select) -> (times, ids, values)

Rows of a simulation table (`species_ids`, then `t_<time>` columns) whose id passes `select`.
The file is memory-mapped, only the selected rows are read. `values[k, j]` is species `k` at
time point `j`.
"""
function read_species(path, select)
    table = Arrow.Table(path)
    rows = findall(select, table.species_ids)
    time_columns = propertynames(table)[2:end]
    times = parse.(Float64, chopprefix.(string.(time_columns), "t_"))
    values = Matrix{Float64}(undef, length(rows), length(time_columns))
    for (j, c) in enumerate(time_columns)
        values[:, j] .= view(getproperty(table, c), rows)
    end
    return times, table.species_ids[rows], values
end

"""
    plot_simulation(tree_configuration, n_node; dt = 0.02, flow_factor = 1.0,
                    zoom = (6.8, 8.0), divider = 7.0, max_terminals = 200, size = (1100, 640))

Concentrations of one saved simulation, like `plot_julia_simulations_with_zoom` (Python):
rows = inflow trees (A, P: marginal edge = what enters the liver), liver tissue (terminal
compartments T_k), outflow trees (V, B: marginal edge = what leaves the liver). Each row has
the whole time span and a zoom of `zoom` (`nothing`: no zoom). At most `max_terminals`
terminal curves are drawn (evenly chosen), so that large trees stay readable. `size`: figure
size in pixels (smaller height for slides with more content).
"""
function plot_simulation(tree_configuration, n_node; dt = 0.02, flow_factor = 1.0,
                         zoom = (6.8, 8.0), divider = 7.0, max_terminals = 200,
                         size = (1100, 640),
                         trees = (inflow = ["A", "P"], outflow = ["V", "B"]))
    graph_dir = joinpath(RESULTS_DIR, tree_configuration, "$(tree_configuration)_$(n_node)")
    flow_part = flow_factor == 1.0 ? "" : "_$(flow_factor)Q"
    path(id) = joinpath(graph_dir, "simulations", "$(id)$(flow_part)_simulations_$(dt)_dt.arrow")

    # terminal compartments: evenly chosen subset
    terminal_ids = Arrow.Table(path("T")).species_ids
    all_terminals = findall(startswith("T_"), terminal_ids)
    step = max(1, cld(length(all_terminals), max_terminals))
    chosen = Set(terminal_ids[all_terminals[1:step:end]])
    times, _, terminal_values = read_species(path("T"), in(chosen))

    marginal(id) = read_species(path(id), ==("Marginal"))[3][1, :]
    names = Dict("A" => "Arterial Tree (A)", "P" => "Portal Tree (P)",
                 "V" => "Venous Tree (V)", "B" => "Biliary Tree (B)")
    rows = [
        ("Inflow Vascular Trees", [(marginal(t), TREE_COLORS[t], names[t]) for t in trees.inflow]),
        ("Liver Tissue", [(terminal_values[k, :], TREE_COLORS["T"], "Terminal nodes (T)")
                          for k in axes(terminal_values, 1)]),
        ("Outflow Vascular Trees", [(marginal(t), TREE_COLORS[t], names[t]) for t in trees.outflow]),
    ]

    shown = length(chosen) < length(all_terminals) ? " ($(length(chosen)) of $(length(all_terminals)) terminals shown)" : ""
    fig = Figure(; size)
    Label(fig[0, 1:2], "$(n_node) terminal nodes, $(dt) min synchronization period$(shown)";
          font = :bold, fontsize = 16, tellwidth = false)
    full_axes = Axis[]
    for (r, (title, curves)) in enumerate(rows)
        ax = Axis(fig[r, 1]; title, SUBTITLE_STYLE..., titlealign = :left,
                  xticklabelsvisible = r == length(rows))
        push!(full_axes, ax)
        isnothing(zoom) || vspan!(ax, zoom...; color = ZOOM_COLOR)
        for (y, color, label) in curves
            lines!(ax, times, y; color, label, linewidth = length(curves) > 10 ? 0.6 : 2)
        end
        axislegend(ax; position = :rt, unique = true, framevisible = false, labelsize = 12)
        isnothing(zoom) && continue

        zax = Axis(fig[r, 2]; backgroundcolor = ZOOM_COLOR, xticklabelsvisible = r == length(rows))
        in_zoom = zoom[1] .<= times .<= zoom[2]
        for (y, color, _) in curves
            lines!(zax, times[in_zoom], y[in_zoom]; color, linewidth = length(curves) > 10 ? 0.6 : 2)
        end
        isnothing(divider) || vlines!(zax, divider; color = "#7209b7", linestyle = :dash, linewidth = 1.5)
        xlims!(zax, zoom...)
    end
    linkaxes!(full_axes...)
    xlims!(full_axes[1], extrema(times)...)
    isnothing(zoom) || colsize!(fig.layout, 2, Relative(0.25))
    Label(fig[length(rows) + 1, 1:2], "Time, [min]"; tellwidth = false, LABEL_STYLE..., fontsize = 16)
    Label(fig[1:length(rows), 0], "Concentration"; rotation = π / 2, tellheight = false, LABEL_STYLE..., fontsize = 16)
    return fig
end
