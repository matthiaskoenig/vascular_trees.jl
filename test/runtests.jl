using Test
using VascularTrees

const DATA_DIR = joinpath(@__DIR__, "data")

# Run the full test suite from the package root:
#   julia --project=. -e 'using Pkg; Pkg.test()'
# or from the Julia REPL in package mode (press `]`), with the package environment active:
#   pkg> test
#
# Test-only dependencies (Aqua, Test, ...) belong in test/Project.toml,
# not in the main Project.toml (otherwise Aqua reports them as stale dependencies).
#
# To run only some of the tests, pass their names as test_args, e.g.
#   julia --project=. -e 'using Pkg; Pkg.test(test_args=["mesh"])'
# Names: definitions, vtkhdf, mesh, closest, pipeline, aqua. No args runs everything.
const SELECTED = ARGS
run_test(name) = isempty(SELECTED) || name in SELECTED

@testset "VascularTrees" begin
    run_test("definitions") && @testset "definitions and options" include("test_definitions.jl")
    run_test("vtkhdf")      && @testset "write_vtkhdf" include("test_vtkhdf.jl")
    run_test("mesh")        && @testset "mesh partition" include("test_mesh_partition.jl")
    run_test("closest")     && @testset "find closest nodes" include("test_closest_nodes.jl")
    run_test("pipeline")    && @testset "pipeline Rectangle_quad_10" include("test_pipeline.jl")
    run_test("aqua")        && @testset "Aqua" include("test_aqua.jl")
end
