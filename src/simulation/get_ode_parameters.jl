"""
    ODEparameters

Load the `.arrow` tables written by the trees processing and turn them into the parameters
(`p`) and initial values (`u0`) of the ODE functions in pharmacokinetic_models.jl.

# Input
- `tree_info`: a `TreeStructure`; `tree_info.GRAPH_DIR` points to the graph folder
  (e.g. `.../Rectangle_quad/Rectangle_quad_10`).
- `graph_subsystem`: id of one part of the system, either a single tree (`"A"`, `"P"`, `"V"`, `"B"`)
  or `"T"` for the terminal nodes. The table read is `<GRAPH_DIR>/graphs/<graph_subsystem>.arrow`.

# Output
- `p`: a `VascularTreeParameters` (single tree) or a `TerminalParameters` (`"T"`),
  see definitions.jl.
- `u0`: zero initial values with the shape of the state of that subsystem.
"""
module ODEparameters

    export get_ODE_parameters,
        get_initial_values

    using DataFrames, Arrow
    using ..Definitions:
        TerminalParameters,
        VascularTreeParameters

    using ..TableUtils:
        load_arrow

    """
        get_ODE_parameters(tree_info, graph_subsystem, flow_scaling_factor) -> p

    Read `<tree_info.GRAPH_DIR>/graphs/<graph_subsystem>.arrow` and build the ODE parameters:
    a `TerminalParameters` for `"T"`, a `VascularTreeParameters` for every other id.
    All flows are multiplied by `flow_scaling_factor` (1.0 = original flow).
    """
    function get_ODE_parameters(tree_info, graph_subsystem::String, flow_scaling_factor::Float64)
        # get path to the graph
        ARROW_PATH::String = normpath(joinpath(tree_info.GRAPH_DIR, "graphs/$(graph_subsystem).arrow"))

        if graph_subsystem != "T"
            p = get_tree_parameters(ARROW_PATH, flow_scaling_factor)
        else
            p = get_terminal_parameters(ARROW_PATH, flow_scaling_factor)
        end

        return p
    end

    """
        get_initial_values(state_template) -> Array{Float64}

    Zero initial values with the same shape as `state_template`: one value per edge for a tree
    (called with `p.ODE_groups`), an `(n_inflow+1) * n_terminals` matrix for the terminal nodes
    (called with `p.flow_values`).
    """
    function get_initial_values(state_template::Array{T}) where {T<:Number}
        u0 = zeros(size(state_template))
        return u0
    end

    """
        get_tree_parameters(ARROW_PATH, flow_scaling_factor) -> VascularTreeParameters

    Parameters of one vessel tree, one entry per edge (row of the table written by
    `build_arrow_df`). `id` and `is_inflow` are taken from the first row; the other columns are
    copied into plain `Vector`s (Arrow columns are read-only views), flows are scaled by
    `flow_scaling_factor`.
    """
    function get_tree_parameters(ARROW_PATH::String, flow_scaling_factor::Float64)
        tree_df = load_arrow(ARROW_PATH)
        p = VascularTreeParameters(
            tree_df.vascular_tree_id[1],
            tree_df.is_inflow[1],
            Vector(tree_df.species_ids),
            Vector(tree_df.flows) .* flow_scaling_factor,
            Vector(tree_df.volumes),
            Vector(tree_df.ODE_groups),
            Vector(Vector.(tree_df.pre_elements)),
            Vector(Vector.(tree_df.post_elements)),
        )

        return p
    end

    """
        get_terminal_parameters(ARROW_PATH, flow_scaling_factor) -> TerminalParameters

    Parameters of the terminal nodes from `T.arrow` (one row per terminal, see
    `assemble_tnodes_table`). The `species_id_*` and `flow_value_*` columns are put back into
    `(n_inflow+1) * n_terminals` blocks: row 1 = terminal / summed outflow, then one row per inflow
    tree. `volumes` is the volume of one terminal (the same for all terminals).
    """
    function get_terminal_parameters(ARROW_PATH::String, flow_scaling_factor::Float64)
        tnodes_df = load_arrow(ARROW_PATH)
        p = TerminalParameters(
            ; id = "T",
            # columns species_id_1..k → rows of a (n_inflow+1) * n_terminals block
            species_ids = permutedims(Matrix{String}(tnodes_df[:, r"^species_id_"])),
            flow_values = permutedims(Matrix{Float64}(tnodes_df[:, r"^flow_value_"])) .* flow_scaling_factor,
            volumes = tnodes_df.volume[1],
        )

        return p
    end

end
