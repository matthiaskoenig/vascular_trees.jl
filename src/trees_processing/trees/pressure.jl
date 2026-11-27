function compute_hydraulic_conductance(d, l)
    # mm⁴ / (Pa·s * mm) -> mm³ / (Pa·s)
    μ = 3.71 * 10^-3 # Pa·s, dynamic viscosity of blood
    c = π * d^4 / (128 * μ * l)
    return c
end

const MMHG_TO_PA = 133.322 # 1 mmHg in Pa (conductances are in mm³/(Pa·s), so pressures must be in Pa)

"""
    calculate_pressure_distribution(tree_info; html_report = false)

Assemble the pressure system of all trees in `tree_info`.
If `html_report = true`, the symbolic matrix is also written as an HTML table to
`<GRAPH_DIR>/reports/pressure_matrix_<graph_id>.html` (open it in a browser).
"""
function calculate_pressure_distribution(tree_info; html_report::Bool = false)
    n_variables = 2
    trees = tree_info.vascular_trees

    tree_dfs = Tuple{String, DataFrame}[] # (tree id, edges df); a Vector keeps load order (Dict order is arbitrary)
    for tree in trees
        ARROW_PATH::String = normpath(joinpath(tree_info.GRAPH_DIR, "graphs/$(tree).arrow"))
        df = load_arrow(ARROW_PATH)
        n_variables += count(!ismissing, df.nodes_ids)
        push!(tree_dfs, (df.vascular_tree_id[1], df))
    end

    right_side = spzeros(n_variables * 2)
    n_rows = n_cols = n_variables * 2 # FIXME: final system size

    # The sparse matrix is collected as a list of entries instead of being written cell by cell
    # (inserting into a SparseMatrixCSC one cell at a time is slow).
    # Entry n means: matrix[rows[n], cols[n]] += values[n]
    rows   = Int[]     # equation index (matrix row) of each entry
    cols   = Int[]     # pressure-unknown index (matrix column) of each entry
    values = Float64[] # conductance coefficient of each entry
    labels = String[]  # same entry as text, e.g. "-c1,2" (only for the printed table)
    pressures = String[] # column names: one pressure unknown per node, e.g. "p3_<tree id>"

    # Columns: first all node pressures (tree by tree), then one inflow Q per inflow tree
    n_pressure_cols = sum(length(df.nodes_ids) for (_, df) in tree_dfs)
    Qinflow_col = Int[]   # matrix column of each inflow tree's Q
    q_labels = String[]   # column names of the Q unknowns, appended after all pressure labels

    # Matrix columns of the pressure p_i at the real node of each marginal edge
    # (inflow: (0, i) -> node i; outflow: (i, 0) -> node i), e.g. for boundary conditions
    inflow_marginal_p_col  = Int[] # one per inflow tree, same order as tree_dfs
    outflow_marginal_p_col = Int[] # one per outflow tree, same order as tree_dfs
    inflow_ids  = String[] # tree id belonging to each entry of inflow_marginal_p_col
    outflow_ids = String[] # tree id belonging to each entry of outflow_marginal_p_col
    # Matrix columns of the leaf-node pressure of every preterminal (group 2) edge, all trees
    # (inflow (i, j) -> leaf j; outflow, reversed, (i, j) -> leaf i)
    preterminal_p_cols = Int[]

    # Each tree gets its own block in the matrix: its equations start after
    # the previous trees' equations, its pressure columns after the previous trees' nodes.
    row_offset = 0
    col_offset = 0
    tree_col_offsets = Int[] # column offset of each tree (same order as tree_dfs), reused after the loop
    for (id, df) in tree_dfs
        push!(tree_col_offsets, col_offset)
        is_inflow = df.is_inflow[1]
        # inflow trees: flow splits into successors; outflow trees: flow merges from predecessors
        neighbors = is_inflow ? df.successors : df.predecessors
        add_equations!(rows, cols, values, labels, df.all_edges, df.conductances, df.ODE_groups,
            neighbors, is_inflow, row_offset, col_offset)

        if is_inflow
            q_col = n_pressure_cols + length(Qinflow_col) + 1 # next free column after all pressures
            # uses THIS tree's offsets: the root edge's row and pressure columns lie in this tree's block
            add_inflow_boundary!(rows, cols, values, labels, df.all_edges, df.conductances, df.ODE_groups,
                neighbors, row_offset, col_offset, q_col)
            push!(Qinflow_col, q_col)
            push!(q_labels, "Q$(id)")
        end

        # column of p_i of the marginal edge, within THIS tree's block
        marginal_p_col = col_offset + marginal_node(df.all_edges, df.ODE_groups, is_inflow)
        push!(is_inflow ? inflow_marginal_p_col : outflow_marginal_p_col, marginal_p_col)
        push!(is_inflow ? inflow_ids : outflow_ids, id)

        # leaf-pressure columns of this tree's preterminal edges, within THIS tree's block
        for (k, group) in enumerate(df.ODE_groups)
            group == 2 || continue # preterminal
            i, j = df.all_edges[k]
            push!(preterminal_p_cols, col_offset + (is_inflow ? j : i))
        end

        n_nodes = length(df.nodes_ids)
        append!(pressures, ("p$(i)_$(id)" for i in 1:n_nodes))
        row_offset += length(df.ODE_groups) # one row per edge
        col_offset += n_nodes               # one column per node
    end
    append!(pressures, q_labels) # label k lands at column n_pressure_cols + k == Qinflow_col[k]

    # Global outflow balance (ONE equation for all outflow trees):
    #   Σ_outflow trees Σ_predecessors (g,i) of the marginal edge  c_gi·(p_g - p_i)  -  Σ_inflow trees Q  =  0
    # All outflow trees write into the same row; `sparse` sums their entries.
    outflow_row = row_offset + 1 # first row after all tree blocks
    for (q_col, q_label) in zip(Qinflow_col, q_labels) # -Q of every inflow tree, pushed once
        push!(rows, outflow_row)
        push!(cols, q_col)
        push!(values, -1.0)
        push!(labels, "-$(q_label)")
    end
    for ((_, df), tree_col_offset) in zip(tree_dfs, tree_col_offsets)
        df.is_inflow[1] && continue
        add_outflow_boundary!(rows, cols, values, labels, df.all_edges, df.conductances, df.ODE_groups,
            df.predecessors, outflow_row, tree_col_offset)
    end
    # Pressure drops: one row per (inflow, outflow) pair, p_in - p_out = (P_in - P_out) from TreePressures
    pressure_drop_row = outflow_row + 1
    terminal_row = add_pressure_drop_boundary!(rows, cols, values, labels, right_side, pressure_drop_row,
        inflow_marginal_p_col, inflow_ids, outflow_marginal_p_col, outflow_ids, TreePressures().pressures_mmHg)

    # Terminal pressures: all leaves of all trees have the same pressure (right side = 0)
    add_terminal_pressure_boundary!(rows, cols, values, labels, terminal_row, preterminal_p_cols)

    # Compact numbering: keep only rows / columns that have at least one entry.
    # Empty rows (edges of groups 2, 3, outflow marginal edges, unused rows of the
    # 2·n_variables FIXME size) and empty columns (unused node numbers) make the matrix singular.
    used_rows = sort!(unique(rows)) # old row numbers that are kept, new row r = used_rows[r]
    used_cols = sort!(unique(cols)) # old column numbers that are kept, new column c = used_cols[c]
    row_map = zeros(Int, n_rows); row_map[used_rows] = 1:length(used_rows) # old row → new row
    col_map = zeros(Int, n_cols); col_map[used_cols] = 1:length(used_cols) # old column → new column
    rows_compact = row_map[rows]
    cols_compact = col_map[cols]
    n_rows_compact = length(used_rows)
    n_cols_compact = length(used_cols)
    unknowns = pressures[used_cols]                    # names of the kept unknowns, e.g. "p3_A", "QA"
    right_side_compact = Vector(right_side)[used_rows] # every row with b ≠ 0 also has matrix entries

    # Build the matrix in one step. Entries with the same (row, col) are summed,
    # e.g. the shared node gets -c_k - Σ c_s from several pushed entries.
    conductances = sparse(rows_compact, cols_compact, values, n_rows_compact, n_cols_compact)
    @info "Pressure system: $(n_rows_compact) equations * $(n_cols_compact) unknowns " *
          "(removed $(n_rows - n_rows_compact) empty rows, $(n_cols - n_cols_compact) empty columns)"

    if html_report
        report_path = normpath(joinpath(tree_info.GRAPH_DIR, "reports", "pressure_matrix_$(tree_info.graph_id).html"))
        write_matrix_html_report(report_path, rows_compact, cols_compact, labels, unknowns,
            n_rows_compact, n_cols_compact)
        @info "Pressure matrix report written to $(report_path)"
    end

    if n_rows_compact != n_cols_compact
        @warn "Pressure system is not square ($(n_rows_compact) * $(n_cols_compact)), not solved."
        return (sol = nothing, unknowns = unknowns)
    end

    sol = solve(LinearProblem(conductances, right_side_compact), UMFPACKFactorization())
    return (sol = sol, unknowns = unknowns) # sol.u[c] is the value of unknowns[c]
end

"""
    marginal_node(edges, groups, is_inflow) -> Int

Index (within the tree) of the real node `i` of the tree's single marginal (group 0) edge:
inflow edges are stored as `(0, i)` → second node, outflow edges as `(i, 0)` → first node.
"""
function marginal_node(edges, groups, is_inflow::Bool)
    k = findfirst(==(0), groups)
    k === nothing && error("No marginal edge (group 0) found in $(is_inflow ? "inflow" : "outflow") tree.")
    a, b = edges[k]
    return is_inflow ? b : a
end

"""
    add_equations!(rows, cols, values, labels, edges, conductances, groups, neighbors, is_inflow, row_offset, col_offset)

Append the flow-conservation equations of one tree as matrix entries
(`rows`, `cols`, `values`) plus their text form (`labels`).

One equation (matrix row) per edge `k = i → j` with conductance `c_k`:

- inflow  (successors `j → h`):   c_k·p_i - c_k·p_j + Σ c_s·(p_h - p_j) = 0
- outflow (predecessors `g → i`): c_k·p_i - c_k·p_j + Σ c_s·(p_i - p_g) = 0

i.e. flow through edge `k` = sum of flows through its neighbor edges `s`.

```
 inflow: flow SPLITS at the shared node j        neighbors = successors of k
                                                 edges[s] = (g, h) with g == j
                            s1
                       +---------> h1
            k          |
     i ----------> j --+                         c_k(p_i - p_j) = c_s1(p_j - p_h1)
                       |    s2                                  + c_s2(p_j - p_h2)
                       +---------> h2


 outflow: flow MERGES at the shared node i       neighbors = predecessors of k
                                                 edges[s] = (g, h) with h == i
            s1
     g1 ---------+
                 |      k
                 +--> i ----------> j            c_k(p_i - p_j) = c_s1(p_g1 - p_i)
                 |                                              + c_s2(p_g2 - p_i)
     g2 ---------+
            s2
```

Matrix row of edge `k` (`row_offset + k`), columns are node pressures (`col_offset + node`):

```
             p_i            p_j            p_far (h or g)
 inflow:    +c_k      -c_k - Σ c_s           +c_s
 outflow:   +c_k + Σ c_s     -c_k            -c_s
```

`row_offset` / `col_offset` shift the tree's block inside the full matrix.
"""
function add_equations!(rows, cols, values, labels, edges, conductances, groups, neighbors,
                        is_inflow::Bool, row_offset::Int, col_offset::Int)
    # Inflow and outflow differ only in signs and in which node is shared with the neighbor edges:
    #                 shared node   far node of neighbor   far coeff   shared coeff
    #   inflow:       j             h (of j → h)           +c_s        -c_s
    #   outflow:      i             g (of g → i)           -c_s        +c_s
    neighbor_sign = is_inflow ? 1 : -1
    far_sign_str, shared_sign_str = is_inflow ? ("+", "-") : ("-", "+")

    for (k, group) in enumerate(groups)
        group == 1 || continue # only "other" edges
        i, j = edges[k]
        c_edge = conductances[k]
        row = row_offset + k

        # own edge term: +c_k·p_i - c_k·p_j
        push!(rows, row, row)
        push!(cols, col_offset + i, col_offset + j)
        push!(values, c_edge, -c_edge)
        push!(labels, "+c$i,$j", "-c$i,$j")

        shared = is_inflow ? j : i # node connecting this edge to its neighbors
        for s in neighbors[k]
            g, h = edges[s]
            if (is_inflow ? g : h) != shared
                @warn "Edge $(i) → $(j) has neighbor $(g) → $(h), which is not a successor/predecessor of the same edge."
            end
            far = is_inflow ? h : g # other end of the neighbor edge
            c_neighbor = conductances[s]

            # neighbor term: ±c_s·p_far ∓ c_s·p_shared
            push!(rows, row, row)
            push!(cols, col_offset + far, col_offset + shared)
            push!(values, neighbor_sign * c_neighbor, -neighbor_sign * c_neighbor)
            push!(labels, far_sign_str * "c$g,$h", shared_sign_str * "c$g,$h")
        end
    end
end

"""
    add_inflow_boundary!(rows, cols, values, labels, edges, conductances, groups, neighbors, row_offset, col_offset, q_col)

Inflow boundary of one inflow tree: the flow `Q` entering through the single marginal
(root, group 0) edge `k = i → j` equals the flow through its successors `j → h`:

    Q = Σ c_s·(p_j - p_h)   ⇔   1·Q + Σ c_s·p_h - Σ c_s·p_j = 0

Written into edge `k`'s row of this tree's block (`row_offset + k`, empty otherwise since
`add_equations!` only fills group-1 edges). `row_offset` / `col_offset` must be THIS tree's
offsets; `q_col` is the matrix column of this tree's `Q`.
"""
function add_inflow_boundary!(rows, cols, values, labels, edges, conductances, groups, neighbors,
                              row_offset::Int, col_offset::Int, q_col::Int)
    for (k, group) in enumerate(groups)
        group == 0 || continue # the single "marginal" (root) edge
        i, j = edges[k]
        row = row_offset + k

        # Q term
        push!(rows, row)
        push!(cols, q_col)
        push!(values, 1.0)
        push!(labels, "+1")

        shared = j # node connecting the root edge to its successors
        for s in neighbors[k]
            g, h = edges[s]
            if g != shared
                @warn "Edge $(i) → $(j) has neighbor $(g) → $(h), which is not a successor of the same edge."
            end
            far = h # other end of the successor edge
            c_neighbor = conductances[s]

            # successor term: +c_s·p_far - c_s·p_shared
            push!(rows, row, row)
            push!(cols, col_offset + far, col_offset + shared)
            push!(values, c_neighbor, -c_neighbor)
            push!(labels, "+c$g,$h", "-c$g,$h")
        end
        return nothing # only one marginal edge per tree
    end
    @warn "No marginal edge (group 0) found in inflow tree."
    return nothing
end

"""
    add_outflow_boundary!(rows, cols, values, labels, edges, conductances, groups, neighbors, row, col_offset)

Add one outflow tree's part of the global outflow balance (the `-Q` terms are added by the caller).
Outflow edges point towards the outlet: the marginal (group 0) edge is `k = i → j` (`j` = outlet,
e.g. `(1, 0)`), its predecessors are `g → i`. The outflow of this tree is written through the
predecessors only — the marginal edge itself (its conductance and node `j`) never enters the matrix:

    Σ_predecessors c_gi·(p_g - p_i)

`row` is the shared outflow row (same for all outflow trees), `col_offset` is THIS tree's column offset.
"""
function add_outflow_boundary!(rows, cols, values, labels, edges, conductances, groups, neighbors,
                               row::Int, col_offset::Int)
    for (k, group) in enumerate(groups)
        group == 0 || continue # the single "marginal" (root) edge
        i, j = edges[k]
        shared = i # node connecting the marginal edge to its predecessors
        for s in neighbors[k]
            g, h = edges[s]
            if h != shared
                @warn "Edge $(i) → $(j) has neighbor $(g) → $(h), which is not a predecessor of the same edge."
            end
            far = g # other end of the predecessor edge
            c_neighbor = conductances[s]

            # predecessor term: +c_s·p_far - c_s·p_shared
            push!(rows, row, row)
            push!(cols, col_offset + far, col_offset + shared)
            push!(values, c_neighbor, -c_neighbor)
            push!(labels, "+c$g,$h", "-c$g,$h")
        end
        return nothing # only one marginal edge per tree
    end
    @warn "No marginal edge (group 0) found in outflow tree."
    return nothing
end

"""
    add_pressure_drop_boundary!(rows, cols, values, labels, right_side, row,
                                inflow_marginal_p_col, inflow_ids, outflow_marginal_p_col, outflow_ids, tree_pressures)

One row per (inflow, outflow) tree pair, e.g. A–V, A–B, P–V, P–B:

    p_in - p_out = (P_in - P_out) · MMHG_TO_PA

`p_in` / `p_out` are the pressure unknowns at the marginal nodes (columns `*_marginal_p_col`),
`P_in` / `P_out` the known tree pressures in mmHg (`tree_pressures`, from `TreePressures`), looked up
by tree id. The known difference goes to `right_side`. Returns the next free row.
"""
function add_pressure_drop_boundary!(rows, cols, values, labels, right_side, row,
                                     inflow_marginal_p_col, inflow_ids, outflow_marginal_p_col, outflow_ids,
                                     tree_pressures)
    for (inflow_col, inflow_id) in zip(inflow_marginal_p_col, inflow_ids)
        for (outflow_col, outflow_id) in zip(outflow_marginal_p_col, outflow_ids)
            # +p_in
            push!(rows, row)
            push!(cols, inflow_col)
            push!(values, 1.0)
            push!(labels, "+1")
            # -p_out
            push!(rows, row)
            push!(cols, outflow_col)
            push!(values, -1.0)
            push!(labels, "-1")
            # known pressure difference of this pair, in Pa
            right_side[row] = (tree_pressures[inflow_id] - tree_pressures[outflow_id]) * MMHG_TO_PA
            row += 1
        end
    end
    return row # next free row
end

"""
    add_terminal_pressure_boundary!(rows, cols, values, labels, row, preterminal_p_cols)

All leaf nodes (ends of preterminal edges) of all trees have one common pressure.
Every leaf is tied to the first one, no matching between trees:

    p_ref - p_m = 0    for every other leaf m      (right side 0, nothing to write)

`preterminal_p_cols` are the matrix columns of the leaf pressures. Returns the next free row.
"""
function add_terminal_pressure_boundary!(rows, cols, values, labels, row, preterminal_p_cols)
    isempty(preterminal_p_cols) && return row
    ref_col = preterminal_p_cols[1]
    for col in @view preterminal_p_cols[2:end]
        push!(rows, row, row)
        push!(cols, ref_col, col)
        push!(values, 1.0, -1.0)
        push!(labels, "+1", "-1")
        row += 1
    end
    return row # next free row
end

# Style of the HTML matrix report:
# small font, tight padding, no line wrapping inside cells,
# column labels (top) and row labels (left) stay visible when scrolling,
# selected rows/columns are tinted blue (box-shadow overlays the cell's own background, e.g. the yellow diagonal),
# the crossing of a selected row and column is tinted darker.
const _MATRIX_REPORT_CSS = raw"""
    body { font-family: sans-serif; }
    table { border-collapse: collapse; }
    table, td, th { font-family: monospace; font-size: 10px; }
    td, th { padding: 1px 3px !important; white-space: nowrap; }
    thead th { position: sticky; top: 0; background: #ffffff; z-index: 1; cursor: pointer; }
    thead th.stubheadLabel { left: 0; z-index: 2; }
    td.rowLabel { position: sticky; left: 0; background: #ffffff; cursor: pointer; }
    .sel-row, .sel-col { box-shadow: inset 0 0 0 9999px rgba(0, 120, 255, 0.15); }
    .sel-row.sel-col   { box-shadow: inset 0 0 0 9999px rgba(0, 120, 255, 0.45); }
    .help { font-size: 12px; color: #555555; margin: 4px 0 8px 0; }
"""

# Click handling for row/column selection in the HTML matrix report.
# PrettyTables' table structure: header <tr class="columnLabelRow"> with <th> per column
# (cell 0 = empty corner), body <tr class="dataRow"> with <td class="rowLabel"> (cell 0) + data cells.
# So a column has the same cellIndex in the header and in every body row.
const _MATRIX_REPORT_JS = raw"""
(function () {
  const table    = document.querySelector("table");
  const header   = table.querySelector("tr.columnLabelRow");
  const dataRows = Array.from(table.querySelectorAll("tbody tr.dataRow"));
  const selCols  = new Set(); // selected column indices (cellIndex)
  const selRows  = new Set(); // selected <tr> elements

  function setCol(c, on) {
    on ? selCols.add(c) : selCols.delete(c);
    header.cells[c].classList.toggle("sel-col", on);
    for (const tr of dataRows) tr.cells[c].classList.toggle("sel-col", on);
  }
  function setRow(tr, on) {
    on ? selRows.add(tr) : selRows.delete(tr);
    for (const td of tr.cells) td.classList.toggle("sel-row", on);
  }
  function clearAll() {
    for (const c of Array.from(selCols)) setCol(c, false);
    for (const tr of Array.from(selRows)) setRow(tr, false);
  }

  // column label: toggle that column
  header.addEventListener("click", e => {
    const th = e.target.closest("th");
    if (th && th.cellIndex > 0) setCol(th.cellIndex, !selCols.has(th.cellIndex));
  });

  // row label: toggle that row; data cell: toggle its row AND column together
  table.tBodies[0].addEventListener("click", e => {
    const td = e.target.closest("td");
    if (!td) return;
    const tr = td.parentElement;
    if (td.cellIndex === 0) {
      setRow(tr, !selRows.has(tr));
    } else {
      const on = !(selRows.has(tr) && selCols.has(td.cellIndex));
      setRow(tr, on);
      setCol(td.cellIndex, on);
    }
  });

  // Esc: clear all selections
  document.addEventListener("keydown", e => { if (e.key === "Escape") clearAll(); });
})();
"""

"""
    write_matrix_html_report(path, rows, cols, labels, pressures, n_rows, n_cols)

Write the symbolic matrix (cell = concatenated `labels` of its entries, `"0"` if empty)
as a stand-alone HTML page to `path`. Compact layout; the diagonal has a yellow background,
zero cells are greyed out, non-zero cells are bold.

Interactive selection in the browser: click a column label to highlight the column,
a row label to highlight the row, a data cell to highlight its row and column; `Esc` clears.
"""
function write_matrix_html_report(path, rows, cols, labels, pressures, n_rows, n_cols)
    # Text version of the matrix; labels of the same cell are concatenated, e.g. "-c1,2-c2,3"
    matrix_str = fill("0", n_rows, n_cols)
    for (row, col, label) in zip(rows, cols, labels)
        cell = matrix_str[row, col]
        matrix_str[row, col] = cell == "0" ? label : cell * label
    end

    mkpath(dirname(path))
    open(path, "w") do io
        # The page is written by hand (instead of PrettyTables' `stand_alone = true`)
        # so that the selection script can be added after the table.
        println(io, "<!DOCTYPE html>\n<html>\n<head>\n<meta charset=\"UTF-8\">")
        println(io, "<style>", _MATRIX_REPORT_CSS, "</style>\n</head>\n<body>")
        println(io, "<div class=\"help\">Click a column label / row label to highlight it, ",
            "a cell to highlight its row and column. Esc clears the selection.</div>")
        pretty_table(io, matrix_str;
            backend       = :html, # writes only the <table> (stand_alone = false)
            title         = "Pressure system: $(n_rows) equations * $(n_cols) unknowns",
            column_labels = vcat(pressures, fill("", n_cols - length(pressures))), # pad to table width
            row_labels    = ["eq$i" for i in 1:n_rows],
            alignment     = :c,
            # Only the FIRST matching highlighter is applied to a cell, so the diagonal goes first
            highlighters  = [
                HtmlHighlighter((d, i, j) -> i == j,  # diagonal: yellow background
                    (h, d, i, j) -> d[i, j] == "0" ?
                        ["background-color" => "#ffe680", "color" => "#999999"] :
                        ["background-color" => "#ffe680", "font-weight" => "bold"]),
                HtmlHighlighter((d, i, j) -> d[i, j] == "0", ["color" => "#cccccc"]),  # grey out zeros
                HtmlHighlighter((d, i, j) -> d[i, j] != "0", ["font-weight" => "bold"]),
            ],
        )
        println(io, "<script>", _MATRIX_REPORT_JS, "</script>\n</body>\n</html>")
    end
    return path
end