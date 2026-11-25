function voxelize_mesh(mesh_path::Vector{String}, n_terminals::Int; check_partial_voxels::Bool = false)
    (; V, vertices, triangles) = compute_mesh_volume(mesh_path)

    volume_per_terminal = V / n_terminals
    average_spacing = (volume_per_terminal)^(1/3)
    voxel_size = average_spacing / 3

    grid = create_grid(vertices, voxel_size)
    inside, columns_with_holes = voxelize(vertices, triangles, grid)

    partly_outside, partly_inside = check_partial_voxels ?
        find_partial_voxels(vertices, triangles, grid, inside) : (nothing, nothing)

    return (
        volume = V,
        inside = inside,
        columns_with_holes = columns_with_holes,
        grid = grid,
        average_spacing = average_spacing,
        partly_outside = partly_outside,
        partly_inside = partly_inside,
    )
end

"""
    find_partial_voxels(X, F, g, inside) -> (partly_outside, partly_inside)

Finds the voxels that are cut by the surface (X, F), using their 8 corners.

`voxelize` decides inside/outside from the voxel center only. Here the corners are
voxelized too (a grid shifted by -h/2 with one more point per axis, so its "centers"
are the corners of the voxels of `g`). Then:
- `partly_outside[i, j, k]`: the voxel is marked inside, but at least one corner is outside;
- `partly_inside[i, j, k]`: the voxel is marked outside, but at least one corner is inside.

Surface features thinner than one voxel can be missed.
"""
function find_partial_voxels(X, F, g::VoxelGrid, inside::BitArray{3})
    corner_grid = VoxelGrid(g.origin .- g.h / 2, g.h, g.dims .+ 1)
    corner_inside, _ = voxelize(X, F, corner_grid)

    nx, ny, nz = g.dims
    partly_outside = falses(nx, ny, nz)
    partly_inside = falses(nx, ny, nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        corners = @view corner_inside[i:i+1, j:j+1, k:k+1]
        if inside[i, j, k]
            partly_outside[i, j, k] = !all(corners)
        else
            partly_inside[i, j, k] = any(corners)
        end
    end
    return partly_outside, partly_inside
end

function compute_mesh_volume(mesh_path::Vector{String})
    # FIXME: check units, switch to manual calculation if necessary
    isempty(mesh_path) && error("No mesh file given")
    length(mesh_path) > 1 && @warn "Multiple meshes given, using the first one: $(mesh_path[1])"
    isfile(mesh_path[1]) || error("Mesh file not found: $(mesh_path[1])")
    mesh = load(mesh_path[1])
    V = volume(mesh)
    pts = coordinates(mesh) # vector of vertices, each a Point3f (x, y, z)
    tri = faces(mesh) # vector of triangles, each a Face3f (i, j, k) with indices into pts
    return (V = V, vertices = Float64[p[d] for d in 1:3, p in pts], triangles = Int[Int(GeometryBasics.value(t[k])) for k in 1:3, t in tri])
end

function create_grid(X, h)
    lo = vec(minimum(X; dims = 2)) .- h
    hi = vec(maximum(X; dims = 2)) .+ h
    dims = ntuple(d -> ceil(Int, (hi[d] - lo[d]) / h), 3)

    return VoxelGrid(Tuple(lo), h, dims)
end

# Tie-break for a column that passes exactly through a triangle edge: of the two
# directions u→v and v→u, exactly one is accepted, so a column through an edge
# shared by two triangles is counted once, not twice or zero times.
owns_edge(dx, dy) = dy < 0 || (dy == 0 && dx > 0)

"""
    voxelize(X, F, g) -> (inside, columns_with_holes)

Marks the voxels whose centers are inside the closed surface (X, F).

For every column of voxels along z, it finds the heights where the column crosses
the surface. Going up a column, the crossings alternate entering and leaving V,
so the voxels between the 1st and 2nd crossing are inside, between the 3rd and
4th, and so on.

`columns_with_holes` counts columns with an odd number of crossings; that happens only if
the mesh has holes. Those columns are left empty.
"""
function voxelize(X, F, g::VoxelGrid)
    nx, ny, nz = g.dims
    ox, oy, oz = g.origin
    h = g.h
    crossings = [Float64[] for _ in 1:nx, _ in 1:ny]

    for t in axes(F, 2)
        ax, ay, az = X[1, F[1, t]], X[2, F[1, t]], X[3, F[1, t]]
        bx, by, bz = X[1, F[2, t]], X[2, F[2, t]], X[3, F[2, t]]
        cx, cy, cz = X[1, F[3, t]], X[2, F[3, t]], X[3, F[3, t]]

        # Twice the signed area of the triangle seen from above. Zero means the
        # triangle is vertical, so no column crosses it.
        area2 = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        area2 == 0 && continue
        if area2 < 0
            bx, by, bz, cx, cy, cz = cx, cy, cz, bx, by, bz
            area2 = -area2
        end

        # Columns whose centers fall inside the triangle's bounding box (seen from above).
        i_lo = max(1, ceil(Int, (min(ax, bx, cx) - ox) / h + 0.5))
        i_hi = min(nx, floor(Int, (max(ax, bx, cx) - ox) / h + 0.5))
        j_lo = max(1, ceil(Int, (min(ay, by, cy) - oy) / h + 0.5))
        j_hi = min(ny, floor(Int, (max(ay, by, cy) - oy) / h + 0.5))

        for j in j_lo:j_hi, i in i_lo:i_hi
            px = ox + (i - 0.5) * h
            py = oy + (j - 0.5) * h

            # Barycentric weights (times area2) of the column center.
            wa = (cx - bx) * (py - by) - (cy - by) * (px - bx)
            wb = (ax - cx) * (py - cy) - (ay - cy) * (px - cx)
            wc = (bx - ax) * (py - ay) - (by - ay) * (px - ax)

            (wa > 0 || (wa == 0 && owns_edge(cx - bx, cy - by))) || continue
            (wb > 0 || (wb == 0 && owns_edge(ax - cx, ay - cy))) || continue
            (wc > 0 || (wc == 0 && owns_edge(bx - ax, by - ay))) || continue

            push!(crossings[i, j], (wa * az + wb * bz + wc * cz) / area2)
        end
    end

    inside = falses(nx, ny, nz)
    columns_with_holes = 0
    for j in 1:ny, i in 1:nx
        zs = crossings[i, j]
        isempty(zs) && continue
        if isodd(length(zs))
            columns_with_holes += 1
            continue
        end
        sort!(zs)
        for q in 1:2:length(zs)
            k_lo = max(1, ceil(Int, (zs[q] - oz) / h + 0.5))
            k_hi = min(nz, floor(Int, (zs[q+1] - oz) / h + 0.5))
            for k in k_lo:k_hi
                inside[i, j, k] = true
            end
        end
    end
    return inside, columns_with_holes
end