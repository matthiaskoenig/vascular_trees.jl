"""
    ReadTreeFiles

Read vascular tree data from CSV files produced by the SyntheticVascularTree generator
and assemble them into DataFrames ready for graph construction.

## Edge attribute convention

In the SyntheticVascularTree format, all attributes are stored on *nodes*, not edges.
The file writer therefore iterates over nodes (nnodes) rather than edges (nedges).
Attributes of a node describe the edge *arriving at* that node (the "target" edge):

    1 → 2 → 3

- Node 1 (root): has `flow`, `radius` — no `length` or `pressure_drop` (no incoming edge).
- Node 2: has `flow`, `radius`, `length`, `pressure_drop` → these values characterise edge (1, 2).

Inflow/outflow direction is irrelevant at this stage; every tree is read as flowing from
the highest-index node to the lowest. Direction reversal is handled downstream.
"""
module ReadTreeFiles
    using CSV, 
        DataFrames, 
        DataFramesMeta, 
        Parameters, 
        Revise

    export read_edges, 
        read_graph_skeleton, 
        read_edge_attributes, 
        read_nodes_attributes

    MM3_PER_L = 1e6
    S_PER_MIN = 60
    """
        read_edges(GRAPH_PATH, EDGES_PATH) -> DataFrame

    Load a complete edge table by joining the graph skeleton (connectivity) with edge
    attributes (physical quantities).

    The two CSVs are joined on `:target_id`. Missing values are not allowed in the
    result and will raise an error if any join keys are unmatched.

    # Arguments
    - `GRAPH_PATH`: path to the CSV containing edge connectivity (source/target node ids).
    - `EDGES_PATH`: path to the CSV containing per-edge physical attributes.

    # Returns
    A `DataFrame` with columns from both sources, with units converted as described in
    [`read_edge_attributes`](@ref).
    """
    function read_edges(GRAPH_PATH::String, EDGES_PATH::String)::DataFrame
        edges_df::DataFrame = read_graph_skeleton(GRAPH_PATH)
        edges_attrib_df::DataFrame = read_edge_attributes(EDGES_PATH)
        leftjoin!(edges_df, edges_attrib_df, on = :target_id)
        disallowmissing!(edges_df)

        return edges_df
    end

    """
        read_graph_skeleton(GRAPH_PATH) -> DataFrame

    Read edge connectivity from a Graphs.jl LGFormat file and return it as a `DataFrame`
    with columns `[:source_id, :target_id]`.

    The first line of an LGFormat file is a metadata header not column names, so it is skipped.

    # Arguments
    - `GRAPH_PATH`: path to the `.lg` file with source and target node ids.
    """
    function read_graph_skeleton(GRAPH_PATH::String)::DataFrame
        graph_skeleton::DataFrame = CSV.read(GRAPH_PATH, DataFrame;
                                             header=[:source_id, :target_id], skipto=2, delim=',',
                                             types=Dict(:source_id => Int, :target_id => Int))

        return graph_skeleton
    end

    """
        read_edge_attributes(EDGES_PATH) -> DataFrame

    Read per-edge physical attributes from CSV, rename columns to short canonical names,
    and convert units to SI-adjacent values used throughout the pipeline:

    | Original column            | Renamed to      | Unit conversion              |
    |:---------------------------|:----------------|:-----------------------------|
    | `edge_idx`                 | `target_id`     | —                            |
    | `leaf_count`               | `leaf`          | —                            |
    | `radius_in_mm`             | `radius`        | mm (unchanged)               |
    | `flow_in_mm3_per_s`        | `flows`         | mm³/s → L/min                |
    | `length_in_mm`             | `length`        | mm (unchanged)               |
    | `pressure_drop_in_kg_per_mm_s2` | `pressure_drop` | kg/(mm·s²) (unchanged)  |

    A new column `volumes` (in L) is computed as `π · radius² · length / 10⁶`.

    # Arguments
    - `EDGES_PATH`: path to the CSV file with edge attributes.
    """
    function read_edge_attributes(EDGES_PATH::String)::DataFrame
        edges_attrib_df::DataFrame = CSV.read(
            EDGES_PATH, 
            DataFrame,
            types=Dict(
                :edge_idx => Int,
                :leaf_count => Int,
                :radius_in_mm => Float64,
                :flow_in_mm3_per_s => Float64,
                :length_in_mm => Float64,
                :pressure_drop_in_kg_per_mm_s2 => Float64
            ); 
            delim=","
        )
        @chain edges_attrib_df begin
            @rename! begin
                :target_id = :edge_idx
                :leaf = :leaf_count
                :radius = :radius_in_mm
                :flows = :flow_in_mm3_per_s
                :length = :length_in_mm
                :pressure_drop = :pressure_drop_in_kg_per_mm_s2
            end
            @transform! begin
                :flows = (:flows ./ MM3_PER_L * S_PER_MIN)           # mm³/s → L/min
                :volumes = π .* :radius .^ 2 .* :length ./ MM3_PER_L  # mm³ → L
            end
        end

        return edges_attrib_df
    end

    """
        read_nodes(NODES_PATH) -> DataFrame

    Read node spatial coordinates from CSV and return a `DataFrame` with columns
    `[:ids, :x, :y, :z]` (coordinates in mm, unchanged from source).

    # Arguments
    - `NODES_PATH`: path to the CSV file with node coordinates.
    """
    function read_nodes(NODES_PATH::String)::DataFrame
        nodes_df::DataFrame = CSV.read(NODES_PATH, DataFrame)
        @chain nodes_df begin
            @rename! begin
                :ids = :node_idx
                :x = :x_coord_in_mm
                :y = :y_coord_in_mm
                :z = :z_coord_in_mm
            end
        end

        return nodes_df
    end

end