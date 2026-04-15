"""
Top-level test runner for tools/.

Run from repo root:
    julia tools/test/runtests.jl
"""

using Test

@testset "ExtXYZWriter" begin
    include("test_ExtXYZ.jl")
end

@testset "VaspParser" begin
    include("test_VaspParser.jl")
end
