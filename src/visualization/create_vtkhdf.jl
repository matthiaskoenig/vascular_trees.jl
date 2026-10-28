module Visualization

    export write_vtkhdf

    using DataFrames,
        HDF5

    using ..TableUtils:
        load_arrow,
        decompose_lists

    function write_vtkhdf(tree_df, simulations_df, VTK_PATH)
        nodes_coord = collect(skipmissing(tree_df.nodes_coordinates))
        mat = Matrix(reinterpret(reshape, Float64, nodes_coord))

        nodes_row_idx = Dict(node_id => row_idx - 1 for (row_idx, node_id) in enumerate(skipmissing(tree_df.nodes_ids)))

        keep = [edge[1] != 0 && edge[2] != 0 for edge in tree_df.all_edges]
        edges = tree_df.all_edges[keep]
        species = tree_df.species_ids[keep]
        n_cells = count(keep)

        connectivity, offset = decompose_lists(edges, 0)
        for (kn, node_idx) in enumerate(connectivity)
            connectivity[kn] = nodes_row_idx[node_idx]
        end

        n_steps = ncol(simulations_df) - 1

        time_values = parse.(Float64, chopprefix.(names(simulations_df)[2:end], "t_"))

        if simulations_df.species_ids != species
            idx = indexin(species, simulations_df.species_ids)
            @assert !any(isnothing, idx) "some species missing in simulation results"
            simulations_df = simulations_df[Int.(idx), :]
        end

        flows = collect(Float64.(tree_df.flows)[keep])
        volumes = collect(Float64.(tree_df.volumes)[keep])

        type_name = "PolyData"
        dt = HDF5.Datatype(HDF5.API.h5t_copy(HDF5.API.H5T_C_S1))
        HDF5.API.h5t_set_size(dt, sizeof(type_name))
        HDF5.API.h5t_set_cset(dt, HDF5.API.H5T_CSET_ASCII)

        h5open(VTK_PATH, "w") do file
            g = create_group(file, "VTKHDF")
            attributes(g)["Version"] = Int32[2, 3]
            attr = create_attribute(g, "Type", dt, HDF5.dataspace(()))
            write_attribute(attr, dt, type_name)

            g["NumberOfPoints"] = fill(length(nodes_row_idx), n_steps)
            g["Points"] = mat

            lines = create_group(g, "Lines")
            lines["NumberOfCells"] = fill(n_cells, n_steps)
            lines["NumberOfConnectivityIds"] = fill(n_cells * 2, n_steps)
            lines["Offsets"] = offset
            lines["Connectivity"] = connectivity

            for add_group_id in ["Vertices", "Polygons", "Strips"]
                add_group = create_group(g, add_group_id)
                add_group["NumberOfCells"] = zeros(Int64, n_steps)
                add_group["NumberOfConnectivityIds"] =  zeros(Int64, n_steps)
                add_group["Offsets"] =  zeros(Int64, n_steps)
                add_group["Connectivity"] = Int64[]
            end

            cell_data = create_group(g, "CellData")
            cell_data["Flow"] = flows
            cell_data["Volume"] = volumes
            cell_data["Concentration"] = vec(Float64.(Matrix(simulations_df[:, 2:end])))

            steps = create_group(g, "Steps")
            attributes(steps)["NSteps"] = n_steps
            steps["Values"] = time_values
            steps["PartOffsets"] = zeros(Int64, n_steps)
            steps["NumberOfParts"] = ones(Int64, n_steps)
            steps["PointOffsets"] = zeros(Int64, n_steps)
            steps["CellOffsets"] = zeros(Int64, 4, n_steps)
            steps["ConnectivityIdOffsets"] = zeros(Int64, 4, n_steps)

            cell_offset = create_group(steps, "CellDataOffsets")
            cell_offset["Flow"] = zeros(Int64, n_steps)
            cell_offset["Volume"] = zeros(Int64, n_steps)
            cell_offset["Concentration"] = Int64.((0:n_steps-1) .* n_cells)
        end
    end

end