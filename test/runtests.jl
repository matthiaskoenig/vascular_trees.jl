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
@testset "VascularTrees" begin
    @testset "definitions and options" include("test_definitions.jl")
    @testset "write_vtkhdf" include("test_vtkhdf.jl")
    @testset "find closest nodes" include("test_closest_nodes.jl")
    @testset "pipeline Rectangle_quad_10" include("test_pipeline.jl")
    @testset "Aqua" include("test_aqua.jl")
end
