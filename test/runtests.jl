using Test
using VascularTrees

const DATA_DIR = joinpath(@__DIR__, "data")

@testset "VascularTrees" begin
    @testset "definitions and options" include("test_definitions.jl")
    @testset "write_vtkhdf" include("test_vtkhdf.jl")
    @testset "pipeline Rectangle_quad_10" include("test_pipeline.jl")
    @testset "Aqua" include("test_aqua.jl")
end
