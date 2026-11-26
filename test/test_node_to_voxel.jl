# Run only this test:
#   julia --project=. -e 'using Pkg; Pkg.test(test_args=["voxel"])'

using Random
using VascularTrees.TreesProcessing: voxelize_mesh, connect_nodes_to_voxels, neighbors

# Mesh: test/data/rectangle.ply, closed box 100 x 100 x 20 (same as in test_mesh_partition.jl).
tree_info_voxel = TreeStructure(;
    tree_configuration = "Rectangle_quad",
    n_node = 10,
    tree_components = TreeDefinitions().vascular_trees["Rectangle_quad"],
    CONFIG_DIR = DATA_DIR,
)

"""Voxelizes the box for N = size(P, 2) points and connects the points (P: 3 x N) to the voxels."""
function run_connection(P::AbstractMatrix{Float64})
    N = size(P, 2)
    v = voxelize_mesh(tree_info_voxel.MESH, N)
    r = connect_nodes_to_voxels(v.grid, P, v.inside)
    return (; r..., N, inside = v.inside, volume = abs(v.volume))
end

# The voxel size is tied to N (h = (V/N)^(1/3) / 3). With few points the 20-thick box is
# only ~3 voxels high, so the voxelized volume can be off by ~10-15 %.
const VOLUME_RTOL_FEW_POINTS = 0.2
const VOLUME_RTOL_MANY_POINTS = 0.1

"""Prints the main numbers of one case (shown even when all tests pass)."""
function report(name, r)
    N = r.N
    n_orphans = length(r.orphans)
    orphan_share = round(n_orphans / N, digits = 2)
    n_outside = length(r.outside)
    n_without_neighbors = count(==(0), diff(r.offsets))
    volume_ratio = round(sum(r.svolume) / r.volume, digits = 3)
    @info(name, N, n_orphans, orphan_share, n_outside, n_without_neighbors, volume_ratio)
end

"""Checks that must hold for any set of points."""
function check_always(r; volume_rtol)
    @test count(==(0), diff(r.offsets)) == 0           # every point has at least one neighbor
    @test sum(r.svolume) / r.volume ≈ 1 rtol = volume_rtol
    @test all((r.owner .== 0) .== .!r.inside)          # owner is 0 exactly outside V
    @test all(0 .<= r.owner .<= r.N)
    for i in 1:r.N                                     # neighbor lists: no self-loops, symmetric
        nb = neighbors(r.offsets, r.nbrs, i)
        @test i ∉ nb
        @test all(i ∈ neighbors(r.offsets, r.nbrs, j) for j in nb)
    end
    @test all(r.svolume[r.orphans] .== 0)
end

# ---- 1. Regular 4 x 4 x 2 lattice: each point is the center of a 25 x 25 x 10 block ----
P_lattice = hcat([[x, y, z] for x in 12.5:25:87.5, y in 12.5:25:87.5, z in (5.0, 15.0)]...)  # 3 x 32
lattice_index = LinearIndices((4, 4, 2))  # (ix, iy, iz) -> column of P_lattice

@testset "regular lattice" begin
    r = run_connection(P_lattice)
    report("regular lattice", r)
    check_always(r; volume_rtol = VOLUME_RTOL_FEW_POINTS)
    @test isempty(r.orphans)
    @test isempty(r.outside)
    @test all(r.svolume .> 0)
    # Per-cell volumes are NOT compared with 6250: with ~3 voxels across z, the upper and
    # lower layers get 1 and 2 voxel layers, so single cells differ a lot from the ideal.

    # Every point is a neighbor of its lattice neighbors along the axes.
    for I in CartesianIndices(lattice_index)
        i = lattice_index[I]
        nb = neighbors(r.offsets, r.nbrs, i)
        for d in (CartesianIndex(1, 0, 0), CartesianIndex(0, 1, 0), CartesianIndex(0, 0, 1))
            J = I + d
            checkbounds(Bool, lattice_index, J) || continue
            @test lattice_index[J] ∈ nb
        end
    end
end

# ---- 2. Random points in the box (fixed seed): realistic, irregular cells ----
@testset "random points" begin
    rng = Xoshiro(42)
    N = 200
    P_random = vcat(100 .* rand(rng, 1, N), 100 .* rand(rng, 1, N), 20 .* rand(rng, 1, N))  # 3 x N
    r = run_connection(P_random)
    report("random points", r)
    check_always(r; volume_rtol = VOLUME_RTOL_MANY_POINTS)
    @test length(r.orphans) / r.N < 0.1
    @test length(r.outside) <= 2   # only points within rounding error of a face could land outside
end

# ---- 3. Edge cases: lattice + points that must trigger `outside` and `orphans` ----
@testset "edge cases" begin
    i_above, i_twin, i_face = 33, 34, 35
    P_edge = hcat(P_lattice,
        [50.0, 50.0, 25.0],    # 33: above the box            -> must be in `outside`
        P_lattice[:, 1],       # 34: exact copy of point 1    -> one of the twins gets no voxels
        [0.0, 50.0, 10.0],     # 35: exactly on the face x = 0 -> inside or just outside, both fine
    )
    r = run_connection(P_edge)
    report("edge cases", r)
    check_always(r; volume_rtol = VOLUME_RTOL_FEW_POINTS)
    @test i_above ∈ r.outside
    @test r.outside ⊆ [i_above, i_face]
    # The twins share all their voxels; the nearest-neighbor search gives them to one of
    # them, the other becomes an orphan and is attached to its twin.
    @test length(intersect(r.orphans, [1, i_twin])) == 1
    orphan = only(intersect(r.orphans, [1, i_twin]))
    twin = orphan == 1 ? i_twin : 1
    @test twin ∈ neighbors(r.offsets, r.nbrs, orphan)
end

# ---- 4. Stress test: hand-written coords, 16 points packed in a 10 x 10 x 10 corner ----
# Most points share a few voxels, so many orphans are expected; only the rules that hold
# for any input are tested, the share of orphans is just reported.
@testset "clustered coords" begin
    P_clustered = hcat(vcat(
        [[0.0, 0.0, 0.0],  [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],  # tree 1: reference tree
        [[0.0, 10.1, 0.0], [0.1, 0.0, 0.0],  [0.0, 0.0, 10.1], [10.1, 0.0, 0.0]],  # tree 2
        [[0.0, 0.0, 9.9],  [9.9, 0.0, 0.0],  [0.0, 0.1, 0.0],  [0.0, 9.9, 0.0]],   # tree 3
        [[10.0, 0.2, 0.0], [0.0, 0.0, 10.2], [0.0, 9.8, 0.0],  [0.0, 0.0, 0.2]],   # tree 4
    )...)  # 3 x 16
    r = run_connection(P_clustered)
    report("clustered coords", r)
    check_always(r; volume_rtol = VOLUME_RTOL_FEW_POINTS)
    @test length(r.outside) <= 2
end
