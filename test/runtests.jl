using HyperRVE
using Gmsh: gmsh
using Test
using Aqua

include("helpers.jl")

@testset "HyperRVE.jl" begin
    @testset "Types (no Gmsh)" begin include("test_types.jl") end
    @testset "Geometry 3D" begin include("test_geometry3d.jl") end
    @testset "Physical groups" begin include("test_groups.jl") end
    @testset "Periodic meshes" begin include("test_periodicity.jl") end
    @testset "Supercells" begin include("test_supercell.jl") end
    @testset "2D cells" begin include("test_2d.jl") end
    @testset "Mesh control" begin include("test_meshcontrol.jl") end
    @testset "Mesh options" begin include("test_meshoptions.jl") end
    @testset "Interface groups" begin include("test_interfaces.jl") end
    @testset "Capsules and fused phases" begin include("test_fuse.jl") end
    @testset "Strut lattices" begin include("test_lattices.jl") end
    @testset "More lattices (BCCZ, FCCZ, 2D)" begin include("test_lattices_more.jl") end
    @testset "Random packings" begin include("test_random.jl") end
    @testset "Polyhedra and polygons" begin include("test_polyhedra.jl") end
    @testset "Voronoi tessellations" begin include("test_voronoi.jl") end
    @testset "Void phases and foams" begin include("test_foams.jl") end
    @testset "Particle chains" begin include("test_chains.jl") end
    @testset "Particle arrays and laminates" begin include("test_arrays.jl") end
    @testset "Auxetic 2D cells" begin include("test_auxetic.jl") end
    @testset "Session" begin include("test_session.jl") end
    @testset "Robustness" begin include("test_robustness.jl") end
    @testset "Aqua" begin Aqua.test_all(HyperRVE) end
end
