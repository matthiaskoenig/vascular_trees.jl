module Benchmarking
    export save_times_as_csv

    using CSV
    using DataFrames
    using TimerOutputs

    using ..Paths: BENCHMARKING_RESULTS_PATH

    """
        timer_row(timer, call_order, iteration) -> NamedTuple

    One row of the running times table from one entry of `TimerOutputs.todict`.
    """
    timer_row(timer, call_order, iteration) = (
        n_calls = timer["n_calls"],
        call_orders = call_order,
        iteration = iteration,
        times_min = timer["time_ns"] / 1e9 / 60,
        allocated_gbytes = timer["allocated_bytes"] / 2^30,
    )

    """
        save_times_as_csv(times, run_info; path = BENCHMARKING_RESULTS_PATH)

    Add the running times measured in `times` to the CSV file `path`.

    Per graph (`run_simulations` times every graph as `<graph_id>` with one inner timer
    `<graph_id>_<k>` per iteration `k`) one row `"Total time (all iterations)"` and one row
    `"Call"` per iteration are written, in the order of the iterations. Every field of
    `run_info` (e.g. solver, tolerances, number of nodes, version) becomes a column with the
    same value in all these rows.

    The columns are matched by name with the rows already in the file: the whole file is read,
    the new rows are added and the file is written again. A column that exists only in the old
    or only in the new rows is empty (`missing`) in the others, so the old rows stay valid when
    the columns change.
    """
    function save_times_as_csv(times::TimerOutput, run_info::NamedTuple;
                               path::String = BENCHMARKING_RESULTS_PATH)
        rows = NamedTuple[]
        for (graph_id, graph_timer) in TimerOutputs.todict(times)["inner_timers"]
            push!(rows, timer_row(graph_timer, "Total time (all iterations)", missing))
            calls = graph_timer["inner_timers"]
            # inner timers are named "<graph_id>_<k>"; a Dict has no order, so sort by k
            iteration(name) = parse(Int, chopprefix(name, "$(graph_id)_"))
            for name in sort!(collect(keys(calls)); by = iteration)
                push!(rows, timer_row(calls[name], "Call", iteration(name)))
            end
        end
        new_rows = DataFrame(rows)
        # scalars become constant columns
        insertcols!(new_rows, pairs(run_info)...)

        table = isfile(path) ? vcat(CSV.read(path, DataFrame), new_rows; cols = :union) : new_rows
        CSV.write(path, table)
    end

end
