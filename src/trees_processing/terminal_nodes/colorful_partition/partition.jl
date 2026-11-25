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

    (; volume, inside, bad_columns, grid, average_spacing, partly_outside, partly_inside) =
        voxelize_mesh(tree_info.MESH, n_terminals; check_partial_voxels)

    if check_partial_voxels
        n_inside = count(inside)
        n_partly_outside = count(partly_outside)
        n_partly_inside = count(partly_inside)
        percent_of_inside = round(100 * n_partly_outside / max(n_inside, 1), digits = 2)
        @info("Voxels cut by the mesh surface",
              n_inside, n_partly_outside, percent_of_inside, n_partly_inside)
    end

    owner = assign_owners(inside, grid, P)

    outside = Int[]
    for i in axes(P, 2)
        v = voxel_of(grid, view(P, :, i))
        if !checkbounds(Bool, inside, v...) || !inside[v...]
            push!(outside, i)
        end
    end

    bad_columns > 0 && @warn "$bad_columns voxel columns crossed the surface an odd number of times; the mesh is probably not closed"
    !isempty(outside) && @warn "$(length(outside)) data points lie outside V (or within one voxel of its surface)"

    return (
        owner = owner,
        grid = grid,
        average_spacing = average_spacing,
        volume = volume,
        outside = outside,
        bad_columns = bad_columns,
        partly_outside = partly_outside,
        partly_inside = partly_inside,
    )
end

find_voxel_center(g::VoxelGrid, i, j, k) =
    (g.origin[1] + (i - 0.5) * g.h,
     g.origin[2] + (j - 0.5) * g.h,
     g.origin[3] + (k - 0.5) * g.h)

"""
    assign_owners(inside, g, P) -> owner

`owner[i, j, k]` is the index of the data point nearest to the center of voxel
(i, j, k), or 0 if the voxel is outside V. Works slice by slice along z, one
slice per thread, so start Julia with several threads (`julia -t auto`).
"""
function assign_owners(inside::BitArray{3}, g::VoxelGrid, P::AbstractMatrix{Float64})
    nx, ny, nz = g.dims
    owner = zeros(Int32, nx, ny, nz)
    tree = KDTree(P)

    Threads.@threads for k in 1:nz
        n_in = count(@view inside[:, :, k])
        n_in == 0 && continue
        centers = Matrix{Float64}(undef, 3, n_in)
        where_ = Vector{CartesianIndex{2}}(undef, n_in)
        q = 0
        for j in 1:ny, i in 1:nx
            inside[i, j, k] || continue
            q += 1
            x, y, z = find_voxel_center(g, i, j, k)
            centers[1, q] = x
            centers[2, q] = y
            centers[3, q] = z
            where_[q] = CartesianIndex(i, j)
        end
        idxs, _ = nn(tree, centers)
        for q in 1:n_in
            owner[where_[q][1], where_[q][2], k] = Int32(idxs[q])
        end
    end
    return owner
end

"""
    voxel_of(g, p) -> (i, j, k)

Voxel that contains the location `p` (may lie outside the grid).
"""
voxel_of(g::VoxelGrid, p) =
    ntuple(d -> floor(Int, (p[d] - g.origin[d]) / g.h) + 1, 3)