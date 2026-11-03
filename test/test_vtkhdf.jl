using DataFrames, HDF5
using VascularTrees.Visualization: write_vtkhdf

# Tiny tree: marginal edge (0, 1), inner edges (1, 2) and (1, 3), terminal self-edge (2, 2).
# The node columns are shorter than the edge columns and padded with `missing`, as in the
# .arrow tables.
graph_df = DataFrame(
    nodes_ids = [1, 2, 3, missing],
    nodes_coordinates = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), missing],
    all_edges = [(0, 1), (1, 2), (1, 3), (2, 2)],
    species_ids = ["X_0_1", "X_1_2", "X_1_3", "X_2_2"],
    flows = [3.0, 2.0, 1.0, 2.0],
    volumes = [0.1, 0.2, 0.3, 0.4],
)

# rows in a different order than the edges, including the marginal edge (not written)
simulations_df = DataFrame(
    "species_ids" => ["X_2_2", "X_0_1", "X_1_3", "X_1_2"],
    "t_0.0" => [3.0, 0.0, 2.0, 1.0],
    "t_1.0" => [30.0, 0.0, 20.0, 10.0],
)

path = joinpath(mktempdir(), "X.vtkhdf")
write_vtkhdf(graph_df, simulations_df, path)

h5open(path, "r") do file
    vtkhdf = file["VTKHDF"]
    @test read_attribute(vtkhdf, "Version") == Int32[2, 3]
    @test read_attribute(vtkhdf, "Type") == "PolyData"

    # geometry: 3 points, 3 lines (marginal edge left out), repeated for both time points
    @test read(vtkhdf["NumberOfPoints"]) == [3, 3]
    @test read(vtkhdf["Points"]) == [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
    @test read(vtkhdf["Lines/NumberOfCells"]) == [3, 3]
    @test read(vtkhdf["Lines/Offsets"]) == [0, 2, 4, 6]
    @test read(vtkhdf["Lines/Connectivity"]) == [0, 1, 0, 2, 1, 1]   # 0-based point indices

    # cell data in the order of the lines, time steps concatenated
    @test read(vtkhdf["CellData/Flow"]) == [2.0, 1.0, 2.0]
    @test read(vtkhdf["CellData/Volume"]) == [0.2, 0.3, 0.4]
    @test read(vtkhdf["CellData/Concentration"]) == [1.0, 2.0, 3.0, 10.0, 20.0, 30.0]

    steps = vtkhdf["Steps"]
    @test read_attribute(steps, "NSteps") == 2
    @test read(steps["Values"]) == [0.0, 1.0]
    @test read(steps["CellDataOffsets/Concentration"]) == [0, 3]
    @test read(steps["CellDataOffsets/Flow"]) == [0, 0]
end

# a species of the tree missing in the simulation results
@test_throws AssertionError write_vtkhdf(
    graph_df,
    simulations_df[simulations_df.species_ids .!= "X_1_3", :],
    joinpath(mktempdir(), "missing.vtkhdf"),
)
