"""
    Visualization

Writes the simulation results of every vascular tree as a time-dependent VTKHDF file
(`PolyData`, VTKHDF version 2.3) that can be opened in ParaView (≥ 5.12).

# Input (per vascular tree, in `<GRAPH_DIR>/`)
- `graphs/<tree>.arrow` - geometry: node ids and coordinates, edges, flows, volumes
  (see `build_arrow_df` in individual_trees.jl).
- `simulations/<tree>[_<f>Q]_simulations_<dt>_dt.arrow` - one row per species, one column per
  stored time point (see `save_solutions` in output.jl).

# Output
- `<GRAPH_DIR>/vtkhdf/<tree>[_<f>Q]_<dt>_dt.vtkhdf` - see [`vtkhdf_path`](@ref).

# Layout of the VTKHDF file
- the geometry is static: `Points` (one per node) and `Lines` (one per edge) are stored once and
  every time step points to them (all geometry offsets in `Steps` are 0);
- `CellData` holds one value per edge: `Flow` and `Volume` (static) and `Concentration`
  (all time steps concatenated, step `k` starts at `(k - 1) * n_lines`);
- the groups `Vertices`, `Polygons` and `Strips` are required by the format, but empty.

The terminal part (`T`) has no geometry of its own and is not written.
"""
module Visualization

    export write_vtkhdf_files,
        write_vtkhdf

    using DataFrames,
        HDF5

    using ..TableUtils:
        load_arrow,
        decompose_lists

    using ..Simulation: simulation_path

    """
        write_vtkhdf_files(tree_info, flow_scaling_factor, dt)

    Write one VTKHDF file for every vascular tree in `tree_info.vascular_trees`, from its
    `.arrow` graph and the simulation with `flow_scaling_factor` and time step `dt`.

    # Arguments
    - `tree_info`: a `TreeStructure` (provides `vascular_trees` and `GRAPH_DIR`).
    - `flow_scaling_factor`, `dt`: select the simulation file, must be the values used in
      run_simulations.jl (`dt = sim_options.dt`).
    """
    function write_vtkhdf_files(tree_info, flow_scaling_factor, dt)
        @info "Writing VTKHDF files..."
        for vascular_tree ∈ tree_info.vascular_trees
            graph_df = load_arrow(joinpath(tree_info.GRAPH_DIR, "graphs", "$(vascular_tree).arrow"))
            simulations_df = load_arrow(
                simulation_path(tree_info, vascular_tree, flow_scaling_factor, dt, "arrow"),
            )
            write_vtkhdf(
                graph_df,
                simulations_df,
                vtkhdf_path(tree_info, vascular_tree, flow_scaling_factor, dt),
            )
        end
    end

    """
        vtkhdf_path(tree_info, vascular_tree, flow_scaling_factor, dt)

    Path of the VTKHDF file of `vascular_tree`, in the `vtkhdf` folder of the graph (created if
    needed). As in [`simulation_path`](@ref), the flow scaling factor is part of the file name
    only if it differs from 1.0.
    """
    function vtkhdf_path(tree_info, vascular_tree, flow_scaling_factor, dt)
        flow_part = flow_scaling_factor == 1.0 ? "" : "_$(flow_scaling_factor)Q"
        VTKHDF_DIR = mkpath(joinpath(tree_info.GRAPH_DIR, "vtkhdf"))
        return joinpath(VTKHDF_DIR, "$(vascular_tree)$(flow_part)_$(dt)_dt.vtkhdf")
    end

    """
        write_vtkhdf(graph_df, simulations_df, VTKHDF_PATH)

    Write one vascular tree as a time-dependent VTKHDF `PolyData` file (layout: see the module
    docstring).

    # Arguments
    - `graph_df`: the `.arrow` table of the tree (`nodes_ids`, `nodes_coordinates`, `all_edges`,
      `species_ids`, `flows`, `volumes`).
    - `simulations_df`: the simulation table of the tree (`species_ids`, then `t_<time>` columns).
    - `VTKHDF_PATH`: output file, overwritten if it exists.

    Edges with node 0 (the marginal edge, whose source is outside the tree) have no point and
    are left out. The rows of `simulations_df` are matched to the edges by `species_ids`.
    """
    function write_vtkhdf(graph_df, simulations_df, VTKHDF_PATH::String)
        # === Points: one per node, 0-based index in the order of nodes_ids ===
        nodes_coordinates = collect(skipmissing(graph_df.nodes_coordinates))
        # 3 x n_points in Julia, stored as (n_points, 3) in HDF5 (column-major vs row-major)
        points = Matrix(reinterpret(reshape, Float64, nodes_coordinates))
        point_idxs = Dict(
            node_id => k - 1 for (k, node_id) ∈ enumerate(skipmissing(graph_df.nodes_ids))
        )
        n_points = length(point_idxs)

        # === Lines: one per edge, without the marginal edge ===
        is_line = [edge[1] != 0 && edge[2] != 0 for edge ∈ graph_df.all_edges]
        n_lines = count(is_line)
        # node ids of the line ends, replaced by point indices below
        connectivity, offsets = decompose_lists(graph_df.all_edges[is_line], 0)
        for (k, node_id) ∈ enumerate(connectivity)
            connectivity[k] = point_idxs[node_id]
        end

        # === Cell data: rows of the simulation in the order of the lines ===
        species_ids = graph_df.species_ids[is_line]
        if simulations_df.species_ids != species_ids
            rows = indexin(species_ids, simulations_df.species_ids)
            @assert !any(isnothing, rows) "some species are missing in the simulation results"
            simulations_df = simulations_df[Int.(rows), :]
        end
        flows = Float64.(graph_df.flows[is_line])
        volumes = Float64.(graph_df.volumes[is_line])
        # column-major: all lines of the first time point, then of the second, ...
        concentrations = vec(Float64.(Matrix(simulations_df[:, 2:end])))

        # === Time points (columns t_<time>) ===
        times = parse.(Float64, chopprefix.(names(simulations_df)[2:end], "t_"))
        n_time_points = length(times)

        h5open(VTKHDF_PATH, "w") do file
            vtkhdf = create_group(file, "VTKHDF")
            attributes(vtkhdf)["Version"] = Int32[2, 3]
            write_ascii_attribute(vtkhdf, "Type", "PolyData")

            # geometry, the counts are repeated for every time point
            vtkhdf["NumberOfPoints"] = fill(n_points, n_time_points)
            vtkhdf["Points"] = points

            lines = create_group(vtkhdf, "Lines")
            lines["NumberOfCells"] = fill(n_lines, n_time_points)
            lines["NumberOfConnectivityIds"] = fill(2 * n_lines, n_time_points)
            lines["Offsets"] = offsets
            lines["Connectivity"] = connectivity

            for group_name ∈ ["Vertices", "Polygons", "Strips"]
                empty_group = create_group(vtkhdf, group_name)
                empty_group["NumberOfCells"] = zeros(Int64, n_time_points)
                empty_group["NumberOfConnectivityIds"] = zeros(Int64, n_time_points)
                empty_group["Offsets"] = zeros(Int64, n_time_points)
                empty_group["Connectivity"] = Int64[]
            end

            cell_data = create_group(vtkhdf, "CellData")
            cell_data["Flow"] = flows
            cell_data["Volume"] = volumes
            cell_data["Concentration"] = concentrations

            # time steps: geometry offsets 0 (static), data offsets per array
            steps = create_group(vtkhdf, "Steps")
            attributes(steps)["NSteps"] = n_time_points
            steps["Values"] = times
            steps["PartOffsets"] = zeros(Int64, n_time_points)
            steps["NumberOfParts"] = ones(Int64, n_time_points)
            steps["PointOffsets"] = zeros(Int64, n_time_points)
            # (n_time_points, 4) in HDF5: Vertices, Lines, Polygons, Strips
            steps["CellOffsets"] = zeros(Int64, 4, n_time_points)
            steps["ConnectivityIdOffsets"] = zeros(Int64, 4, n_time_points)

            cell_data_offsets = create_group(steps, "CellDataOffsets")
            cell_data_offsets["Flow"] = zeros(Int64, n_time_points)
            cell_data_offsets["Volume"] = zeros(Int64, n_time_points)
            cell_data_offsets["Concentration"] = Int64.((0:n_time_points-1) .* n_lines)
        end
    end

    """
        write_ascii_attribute(parent, name, value)

    Write `value` as a scalar fixed-length ASCII string attribute. HDF5.jl writes strings as UTF-8,
    which the VTKHDF reader of ParaView does not accept for the `Type` attribute.
    """
    function write_ascii_attribute(parent, name::String, value::String)
        datatype = HDF5.Datatype(HDF5.API.h5t_copy(HDF5.API.H5T_C_S1))
        HDF5.API.h5t_set_size(datatype, sizeof(value))
        HDF5.API.h5t_set_cset(datatype, HDF5.API.H5T_CSET_ASCII)
        attribute = create_attribute(parent, name, datatype, HDF5.dataspace(()))
        write_attribute(attribute, datatype, value)
    end

end
