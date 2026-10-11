"""
    TreesProcessingHelpers

Helpers for preparing a vascular tree edge table for ODE model construction.

The typical workflow is:

1. [`label_special_edges!`](@ref): add boolean columns `:preterminal`, `:start` and
   `:terminal` that classify existing edges by their position in the tree.
2. [`create_special_edges!`](@ref): append artificial edges that the model needs but the
   tree files do not contain: one self-loop per leaf node (terminal compartment) and one
   marginal input edge feeding the root.

Step 1 must run before step 2, because the edges created in step 2 are found using the
labels from step 1.
"""
module TreesProcessingHelpers

    using CSV,
        DataFrames,
        DataFramesMeta,
        InteractiveUtils,
        Parameters,
        Revise,
        Arrow

    export label_special_edges!,
        create_special_edges!

    using ..TableUtils: selection_from_df

    # Total volume of the terminal (tissue) compartment, a cube with 0.1 cm sides: [cm^3] -> [L].
    # It is split equally between all terminal edges.
    volume_geometry = (0.100 * 0.100 * 0.10) / 1000

    #=================================================================================================================================#
    """
        label_special_edges!(graph_structure::DataFrame)

    Add the boolean columns `:preterminal`, `:start` and `:terminal` to `graph_structure`
    in place. See [`label_preterminal_edges!`](@ref), [`label_start_edges!`](@ref) and
    [`label_terminal_edges!`](@ref) for what each label means.

    `graph_structure` must have the columns `:source_id` and `:target_id`.
    """
    function label_special_edges!(graph_structure::DataFrame)
        label_preterminal_edges!(graph_structure)
        label_start_edges!(graph_structure)
        label_terminal_edges!(graph_structure)
    end

    """
        label_preterminal_edges!(graph_structure::DataFrame)

    Mark edges that end in a leaf node, i.e. edges whose `target_id` never appears as a
    `source_id`. They are called *preterminal* because each one will later be followed by
    a terminal self-loop edge (see [`create_terminal_edges!`](@ref)).
    """
    function label_preterminal_edges!(graph_structure::DataFrame)
        graph_structure[!, :preterminal] =
            .!in.(graph_structure.target_id, [Set(graph_structure.source_id)])
    end

    """
        label_start_edges!(graph_structure::DataFrame)

    Mark edges that leave the root node, i.e. edges whose `source_id` never appears as a
    `target_id`. A valid tree has exactly one such edge.
    """
    function label_start_edges!(graph_structure::DataFrame)
        graph_structure[!, :start] =
            .!in.(graph_structure.source_id, [Set(graph_structure.target_id)])
    end

    """
        label_terminal_edges!(graph_structure::DataFrame)

    Add the `:terminal` column, set to `false` for all existing edges. Only the self-loop
    edges added by [`create_terminal_edges!`](@ref) are terminal.
    """
    function label_terminal_edges!(graph_structure::DataFrame)
        graph_structure[!, :terminal] = [false for _ ∈ 1:nrow(graph_structure)]
    end

    #=================================================================================================================================#
    """
        create_special_edges!(graph_structure)

    Append the artificial edges needed by the ODE model to `graph_structure` in place:
    the terminal self-loops ([`create_terminal_edges!`](@ref)) and the marginal input edge
    ([`create_marginal_edge!`](@ref)).

    Requires the labels from [`label_special_edges!`](@ref) and the id columns
    `:species_ids`, `:flow_ids` and `:volume_ids`.
    """
    function create_special_edges!(graph_structure)
        create_terminal_edges!(graph_structure)
        create_marginal_edge!(graph_structure)
    end

    """
        create_terminal_edges!(graph_structure)

    Append one self-loop edge `(n, n)` for every leaf node `n`, i.e. the target of every
    preterminal edge. These edges represent the terminal (tissue) compartments; their
    attributes are built by [`collect_terminal_edge_attributes`](@ref).
    """
    function create_terminal_edges!(graph_structure)
        tnodes_ids_df = selection_from_df(
            graph_structure,
            (graph_structure.preterminal .== true, :target_id),
        )
        tnodes_attrib_df = collect_terminal_edge_attributes(tnodes_ids_df)
        append!(graph_structure, tnodes_attrib_df)
    end

    """
        create_marginal_edge!(graph_structure::DataFrame)

    Append the marginal edge `(0, root)`, which feeds the input into the tree. Node `0` is
    a virtual node outside the tree. All physical attributes of the edge are zero, except
    `:volumes`, which is copied from the start edge. Its `:species_ids` is `"Marginal"`.

    The row is pushed by position, so its values must stay in this column order:
    `source_id, target_id, leaf, radius, flows, length, pressure_drop, volumes,
    species_ids, flow_ids, volume_ids, preterminal, start, terminal`.
    If columns are added, removed or reordered, this function must be updated too.
    """
    function create_marginal_edge!(graph_structure::DataFrame)
        start_node_id =
            (selection_from_df(graph_structure, (graph_structure.start .== true, :source_id)))[1]
        push!(
            graph_structure,
            [
                0,
                start_node_id,
                0,
                0.0,
                0.0,
                0.0,
                0.0,
                graph_structure[graph_structure.start.==true, :volumes][1],
                "Marginal",
                "",
                "",
                false,
                false,
                false,
            ],
        )
    end

    """
        collect_terminal_edge_attributes(terminal_node_ids::SubArray) -> DataFrame

    Build the rows for the terminal self-loop edges, one per id in `terminal_node_ids`.

    Each terminal edge `k` (numbered from 1) gets:
    - `source_id = target_id =` the leaf node id;
    - the ids `T_k` (species), `QT_k` (flow) and `VT_k` (volume);
    - a volume of `volume_geometry / n_terminals`, so the terminal compartments share the
      total terminal volume equally;
    - zero for `:leaf`, `:radius`, `:flows`, `:length` and `:pressure_drop`;
    - `:terminal = true`, `:preterminal = :start = false`.
    """
    function collect_terminal_edge_attributes(terminal_node_ids::SubArray)::DataFrame
        n_terminals = length(terminal_node_ids)
        species_ids = ["T_$n_terminal" for n_terminal = 1:n_terminals]
        flow_ids = ["QT_$n_terminal" for n_terminal = 1:n_terminals]
        volume_ids = ["VT_$n_terminal" for n_terminal = 1:n_terminals]
        volume_terminal = volume_geometry / n_terminals
        return DataFrame(
            :source_id => terminal_node_ids,
            :target_id => terminal_node_ids,
            :leaf .=> 0,
            ([:radius, :flows, :length, :pressure_drop] .=> 0.0)...,
            :volumes .=> volume_terminal,
            :species_ids => species_ids,
            :flow_ids => flow_ids,
            :volume_ids => volume_ids,
            ([:preterminal, :start] .=> false)...,
            :terminal .=> true,
        )
    end

end
