# The 6 face neighbors of a voxel (no edges or corners).
const NB6 = (CartesianIndex(1, 0, 0), CartesianIndex(-1, 0, 0),
             CartesianIndex(0, 1, 0), CartesianIndex(0, -1, 0),
             CartesianIndex(0, 0, 1), CartesianIndex(0, 0, -1))

"""
    colorful_partition(tree_info, n_terminals, points; check_partial_voxels = false)

Voxelizes the mesh of `tree_info` and assigns every inside voxel to its nearest point.

With `check_partial_voxels = true`, the voxels cut by the mesh surface are also found
(see [`find_partial_voxels`](@ref)), their counts are logged, and they are returned as
`partly_outside` / `partly_inside` (otherwise both are `nothing`).
"""
function colorful_partition(
    tree_info::TreeStructure,
    n_terminals::Int,
    points::Vector{NTuple{3,Float64}};
    check_partial_voxels::Bool = false,
)
    P = [p[d] for d in 1:3, p in points]

    (; volume, inside, columns_with_holes, grid, average_spacing, partly_outside, partly_inside) =
        voxelize_mesh(tree_info.MESH, n_terminals; check_partial_voxels)

    if check_partial_voxels
        n_inside = count(inside)
        n_partly_outside = count(partly_outside)
        n_partly_inside = count(partly_inside)
        percent_of_inside = round(100 * n_partly_outside / max(n_inside, 1), digits = 2)
        @info("Voxels cut by the mesh surface",
              n_inside, n_partly_outside, percent_of_inside, n_partly_inside)
    end

    columns_with_holes > 0 && @warn "$columns_with_holes voxel columns crossed the surface an odd number of times; the mesh is probably not closed"

    (; owner, outside) = connect_nodes_to_voxels(grid, P, inside)

    return (
        owner = owner,
        grid = grid,
        average_spacing = average_spacing,
        volume = volume,
        outside = outside,
        bad_columns = columns_with_holes,
        partly_outside = partly_outside,
        partly_inside = partly_inside,
    )
end
