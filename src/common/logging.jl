"""
    LogUtils

Log formatting of VascularTrees: color-blind friendly colors of the log prefixes
(`log_metafmt`, set as the global logger by `setup_logging` in the run scripts) and logging of
the stages of a run with their running times (`stage`).
"""
module LogUtils
    export log_metafmt, 
        setup_logging, 
        stage

    using StyledStrings
    # Logging is not a dependency of the package; its log levels, macros and loggers are in Base
    using Base.CoreLogging: 
        @logmsg, 
        Debug, 
        Info, 
        Warn, 
        Error, 
        ConsoleLogger, 
        global_logger

    # Okabe-Ito orange and sky blue
    const STAGE_START_COLOR = "#E69F00"
    const STAGE_END_COLOR = "#56B4E9"

    """
        log_metafmt(level, _module, group, id, file, line)

    Prefix of the log messages: color-blind friendly colors (256-color codes closest to the
    Okabe-Ito palette; the level is also written as text) and the source location
    (`@ module file:line`) only for warnings and errors.
    """
    function log_metafmt(level, _module, group, id, file, line)
        color = level < Info ? 245 :   # Debug: grey
                level < Warn ? 74  :   # Info: sky blue
                level < Error ? 178 :  # Warning: orange
                166                    # Error: vermillion
        prefix = (level == Warn ? "Warning" : string(level)) * ":"
        suffix = level < Warn ? "" :
                "@ $(_module) $(file === nothing ? "" : Base.contractuser(file)):$(line)"
        return color, prefix, suffix
    end

    """
        setup_logging(; debug = false) -> previous logger

    Set the global logger to a `ConsoleLogger` with [`log_metafmt`](@ref); `debug = true` also
    shows the Debug-level stages of VascularTrees (see [`stage`](@ref)).

    Called explicitly in the run scripts, not when the package is loaded: it replaces the logger
    of the whole session (REPL, VS Code, tests).
    """
    function setup_logging(; debug::Bool = false)
        debug && (ENV["JULIA_DEBUG"] = "VascularTrees")
        return global_logger(ConsoleLogger(stderr; meta_formatter = log_metafmt))
    end

    """
        stage(f, name; debug = false)

    Log the start of the stage `name`, run `f()` and log its end with the running time, the
    allocated memory and the GC time. Returns the result of `f()`.

    The start is logged before `f()` runs, so a long stage is visible while it is running.
    Main stages are logged at `Info` level; `debug = true` logs at `Debug` level (fine-grained
    stages, e.g. per subsystem), shown only with `ENV["JULIA_DEBUG"] = "VascularTrees"`.

    Colors (Okabe-Ito palette, color-blind friendly; start and end also differ by the symbol):
    start `▶` orange, end `✔` sky blue. They are shown only in a terminal, log files get
    plain text.
    """
    function stage(f, name; debug::Bool = false)
        level = debug ? Debug : Info
        @logmsg level styled"{(fg=$STAGE_START_COLOR):▶} {bold:$name}"
        stats = @timed f()
        @logmsg level styled"$(round(stats.time; digits = 2)) s, \
                             $(Base.format_bytes(stats.bytes)) allocated, \
                             $(round(stats.gctime; digits = 2)) s GC"
        return stats.value
    end

end
