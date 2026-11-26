using RVEs
using Gmsh: gmsh
using Test
using Aqua

include("helpers.jl")

@testset "RVEs.jl" begin
    @testset "Types (no Gmsh)" begin include("test_types.jl") end
    @testset "Geometry 3D" begin include("test_geometry3d.jl") end
    @testset "Physical groups" begin include("test_groups.jl") end
    @testset "Periodic meshes" begin include("test_periodicity.jl") end
    @testset "Session" begin include("test_session.jl") end
    @testset "Robustness" begin include("test_robustness.jl") end
    @testset "Aqua" begin Aqua.test_all(RVEs) end
end
