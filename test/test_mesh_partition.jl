# Run only this test:
#   julia --project=. -e 'using Pkg; Pkg.test(test_args=["mesh"])'
# Run it and also save a picture of the mesh and the voxel grid to test/output/mesh_partition.png:
#   VT_PLOT=1 julia --project=. -e 'using Pkg; Pkg.test(test_args=["mesh"])'

using VascularTrees.TreesProcessing: voxelize_mesh, compute_mesh_volume, find_voxel_center

GRAPH_DIR = joinpath(mktempdir(), "Rectangle_quad_10")
cp(joinpath(DATA_DIR, "Rectangle_quad_10"), GRAPH_DIR)

# CONFIG_DIR points to test/data, so MESH is picked up from there
# (test/data/rectangle.ply: closed box 100 x 100 x 20, volume 200_000).
tree_info = TreeStructure(;
    tree_configuration = "Rectangle_quad",
    n_node = 10,
    tree_components = TreeDefinitions().vascular_trees["Rectangle_quad"],
    CONFIG_DIR = DATA_DIR,
    GRAPH_DIR = GRAPH_DIR,
)

@test tree_info.MESH == [joinpath(DATA_DIR, "rectangle.ply")]

box_volume = 100.0 * 100.0 * 20.0
n_terminals = 1000  # enough voxels across the 20-thick box for the voxel volume to be accurate

result = voxelize_mesh(tree_info.MESH, n_terminals; check_partial_voxels = true)

@test abs(result.volume) ≈ box_volume rtol = 1e-4
@test result.columns_with_holes == 0  # the box is closed
@test any(result.inside)

# The voxels marked as inside should add up to roughly the box volume.
voxel_volume = count(result.inside) * result.grid.h^3
@test voxel_volume ≈ box_volume rtol = 0.05

# Voxels cut by the surface: they exist, but only in the surface layer, and the two kinds don't overlap.
@test any(result.partly_outside)
@test count(result.partly_outside) < count(result.inside)
@test !any(result.partly_outside .& result.partly_inside)
@test !any(result.partly_outside .& .!result.inside)  # partly outside ⊆ inside
@test !any(result.partly_inside .& result.inside)     # partly inside ⊆ outside

# Without the switch, nothing is computed.
@test voxelize_mesh(tree_info.MESH, n_terminals).partly_outside === nothing

# ---- Optional visual check (only with VT_PLOT=1, see the command at the top) ----
# Kept as separate top-level statements so CairoMakie is fully loaded before it is used.
const MAKE_PLOT = get(ENV, "VT_PLOT", "0") == "1"
MAKE_PLOT && @eval using CairoMakie

if MAKE_PLOT
    (; vertices, triangles) = compute_mesh_volume(tree_info.MESH)
    grid = result.grid
    (; inside, partly_outside, partly_inside) = result
    nx, ny, nz = grid.dims
    voxel_centers(mask) = [Point3f(find_voxel_center(grid, I[1], I[2], I[3])) for I in findall(mask)]

    fig = Figure(size = (1400, 800))

    # Left: mesh (transparent), voxel centers, and the grid's bounding box.
    #   orange = inside, red = inside but partly outside, purple = outside but partly inside
    ax3d = Axis3(fig[1, 1]; title = "mesh + voxel centers (h = $(round(grid.h, digits = 3)))",
                 aspect = :data)
    mesh!(ax3d, permutedims(vertices), permutedims(triangles); color = (:steelblue, 0.25))
    scatter!(ax3d, voxel_centers(inside .& .!partly_outside); markersize = 3, color = :orange,
             label = "inside")
    scatter!(ax3d, voxel_centers(partly_outside); markersize = 4, color = :red,
             label = "inside, partly outside")
    scatter!(ax3d, voxel_centers(partly_inside); markersize = 4, color = :purple,
             label = "outside, partly inside")
    grid_lo = grid.origin
    grid_size = grid.h .* grid.dims
    wireframe!(ax3d, Rect3f(Vec3f(grid_lo...), Vec3f(grid_size...)); color = :gray40)
    Legend(fig[1, 3], ax3d)

    # Right: one horizontal slice through the middle of the grid.
    #   0 = outside, 1 = inside, 2 = inside but partly outside, 3 = outside but partly inside
    k_mid = cld(nz, 2)
    z_mid = find_voxel_center(grid, 1, 1, k_mid)[3]
    x_centers = [find_voxel_center(grid, i, 1, 1)[1] for i in 1:nx]
    y_centers = [find_voxel_center(grid, 1, j, 1)[2] for j in 1:ny]
    slice_class = Float32.(inside[:, :, k_mid])
    slice_class[partly_outside[:, :, k_mid]] .= 2
    slice_class[partly_inside[:, :, k_mid]] .= 3
    ax2d = Axis(fig[1, 2]; title = "slice k = $k_mid (z = $(round(z_mid, digits = 2)))",
                xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(ax2d, x_centers, y_centers, slice_class;
             colormap = [:white, :orange, :red, :purple], colorrange = (0, 3))

    # Bottom: statistics.
    n_inside = count(inside)
    n_partly_outside = count(partly_outside)
    n_partly_inside = count(partly_inside)
    mesh_volume = abs(result.volume)
    stats = """
        inside voxels: $n_inside    grid: $(nx) × $(ny) × $(nz)    h = $(round(grid.h, digits = 4))
        inside, partly outside (red): $n_partly_outside ($(round(100 * n_partly_outside / max(n_inside, 1), digits = 2)) % of inside)
        outside, partly inside (purple): $n_partly_inside ($(round(100 * n_partly_inside / max(n_inside, 1), digits = 2)) % of inside)
        voxel volume: $(round(voxel_volume, digits = 1))    mesh volume: $(round(mesh_volume, digits = 1))    rel. error: $(round(100 * (voxel_volume - mesh_volume) / mesh_volume, digits = 2)) %
        bad columns: $(result.columns_with_holes)"""
    Label(fig[2, 1:3], stats; justification = :left, halign = :left, tellwidth = false)

    out_dir = joinpath(@__DIR__, "output")
    mkpath(out_dir)
    out_file = joinpath(out_dir, "mesh_partition.png")
    save(out_file, fig)
    @info "Mesh partition plot saved to $out_file"
end
