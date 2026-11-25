using RVEs
using Gmsh: gmsh
using Test
using Aqua

# Invariant-based tests: exact node counts/coordinates depend on the Gmsh version,
# so we check geometry (phase volumes), topology (physical groups), periodicity of the
# boundary mesh and basic mesh sanity instead.

const OUTDIR = mktempdir()

"""Total OCC volume of the 3D entities in physical group `tag`."""
function group_volume(tag)
    ents = gmsh.model.getEntitiesForPhysicalGroup(3, tag)
    return sum(gmsh.model.occ.getMass(3, e) for e in ents; init=0.0)
end

"""Tag of the physical group of dimension `dim` called `name` (error if missing)."""
function group_tag(dim, name)
    for (d, t) in gmsh.model.getPhysicalGroups(dim)
        gmsh.model.getPhysicalName(d, t) == name && return t
    end
    error("physical group $name (dim $dim) not found")
end

"""Check that every boundary node on `face_min` has a partner on `face_max` shifted by L along `dir`."""
function periodic_nodes_match(face_min, face_max, dir, L; tol=1e-8 * L)
    _, cmin = gmsh.model.mesh.getNodesForPhysicalGroup(2, group_tag(2, face_min))
    _, cmax = gmsh.model.mesh.getNodesForPhysicalGroup(2, group_tag(2, face_max))
    pmin = reshape(cmin, 3, :)
    pmax = reshape(cmax, 3, :)
    size(pmin, 2) > 0 || return false  # never pass vacuously
    size(pmin, 2) == size(pmax, 2) || return false
    for i in axes(pmin, 2)
        p = pmin[:, i]
        p[dir] += L
        any(j -> maximum(abs.(pmax[:, j] .- p)) < tol, axes(pmax, 2)) || return false
    end
    return true
end

function check_mesh_sanity()
    etypes, etags, _ = gmsh.model.mesh.getElements(3)
    @test !isempty(etypes)
    @test sum(length, etags) > 0
    # every tetrahedron has positive volume (minimum scaled Jacobian > 0)
    for tags in etags
        q = gmsh.model.mesh.getElementQualities(tags, "minSJ")
        @test minimum(q) > 0
    end
end

function run_case(check, name, rve, inclusions)
    startGmsh()
    ShowInfo(0)
    try
        model, tagvol = createGmshModel(rve, inclusions, name)
        createMesh!(model, rve, inclusions, tagvol)
        saveMesh(joinpath(OUTDIR, name * ".msh"))
        @test isfile(joinpath(OUTDIR, name * ".msh"))
        check_mesh_sanity()
        # phases fill the cell exactly
        vtot = prod(rve.size)
        @test group_volume(1) + group_volume(2) ≈ vtot rtol = 1e-6
        names = [gmsh.model.getPhysicalName(d, t) for (d, t) in gmsh.model.getPhysicalGroups()]
        @test length(gmsh.model.getEntitiesForPhysicalGroup(0, 7)) == 8  # Corners
        for n in ("Phase0", "Inclusions", "Surfxmin", "Surfxmax", "Surfymin", "Surfymax",
                  "Surfzmin", "Surfzmax", "Corners")
            @test n in names
        end
        if all(==(1), rve.periodicityFlags)
            for (d, x) in enumerate(("x", "y", "z"))
                @test periodic_nodes_match("Surf$(x)min", "Surf$(x)max", d, rve.size[d])
            end
        end
        check()
    finally
        stopGmsh()
    end
end

sphere_vol(r) = 4π * r^3 / 3

@testset "RVEs.jl" begin

    @testset "Simple cubic, non-periodic" begin
        rve = RVE([1.0, 1.0, 1.0], [0, 0, 0], [0.0, 0.0, 0.0], 0.1)
        run_case("Test1", rve, (Sphere([0.5, 0.5, 0.5], 0.2, 2),)) do
            @test group_volume(2) ≈ sphere_vol(0.2) rtol = 1e-6
        end
    end

    @testset "Simple cubic, periodic (single inclusion, not a tuple)" begin
        rve = RVE([1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
        run_case("Test2", rve, Sphere([0.5, 0.5, 0.5], 0.2, 2)) do
            @test group_volume(2) ≈ sphere_vol(0.2) rtol = 1e-6
        end
    end

    @testset "Sphere across a corner + rotated ellipsoid" begin
        rve = RVE([1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
        incs = (Sphere([1.0, 1.0, 0.5], 0.1, 2),
                Ellipsoid([0.5, 0.5, 0.5], [0.25, 0.1, 0.1], [pi / 4, pi / 4], 2))
        run_case("Test3", rve, incs) do
            # corner sphere: 4 periodic quarters = one full sphere
            @test group_volume(3) ≈ sphere_vol(0.1) rtol = 1e-6
            # the ellipsoid is a dilated sphere (B-spline in OCC): volume only accurate to ~1e-4
            @test group_volume(4) ≈ 4π * 0.25 * 0.1 * 0.1 / 3 rtol = 1e-3
        end
    end

    @testset "Body-centred cubic" begin
        rve = RVE([1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
        incs = (Sphere([0.5, 0.5, 0.5], 0.25, 2), Sphere([1.0, 1.0, 1.0], 0.25, 2))
        run_case("Test4", rve, incs) do
            @test group_volume(2) ≈ 2 * sphere_vol(0.25) rtol = 1e-6
        end
    end

    @testset "Face-centred cubic" begin
        rve = RVE([1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
        incs = (Sphere([0.0, 0.5, 0.5], 0.25, 2), Sphere([0.5, 0.0, 0.5], 0.25, 2),
                Sphere([0.5, 0.5, 0.0], 0.25, 2), Sphere([1.0, 1.0, 1.0], 0.25, 2))
        run_case("Test5", rve, incs) do
            @test group_volume(2) ≈ 4 * sphere_vol(0.25) rtol = 1e-6
        end
    end

    @testset "Fibres" begin
        rve = RVE([1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
        incs = (Cylinder([-0.1, 0.5, 0.5], [2.0, 0.0, 0.0], 0.25, 2),
                Cylinder([-0.1, 1.0, 1.0], [2.0, 0.0, 0.0], 0.25, 2))
        run_case("Test6", rve, incs) do
            @test group_volume(2) ≈ 2 * π * 0.25^2 rtol = 1e-6
        end
    end

    @testset "Parametric fibre array" begin
        rve = RVE([1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
        numfib = [4, 4]
        step = rve.size[2:3] / (numfib[1] - 1)
        radius = minimum(step / 2.0 - 0.1 * step)
        incs = Tuple(Cylinder([-0.1, (iy - 1) * step[1], (iz - 1) * step[2]], [2.0, 0.0, 0.0], radius, 2)
                     for iz in 1:numfib[2] for iy in 1:numfib[1])
        run_case("Test7", rve, incs) do
            # 3x3 full fibres per unit cell (16 lattice points, edges shared periodically)
            @test group_volume(2) ≈ 9 * π * radius^2 rtol = 1e-6
        end
    end

    @testset "Inclusion from boolean operations" begin
        rve = RVE([1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
        geo4 = Fuse([Cylinder([-0.1, 0.5, 0.5], [2.0, 0.0, 0.0], 0.1, 2),
                     Cylinder([0.5, -0.1, 0.5], [0.0, 2.0, 0.0], 0.1, 2),
                     Cylinder([0.5, 0.5, -0.1], [0.0, 0.0, 2.0], 0.1, 2)], 2)
        L = 0.38
        geo7 = Intersect(Sphere([0.5, 0.5, 0.5], 0.25, 2),
                         Box([L, L, L], [0.5 - L / 2, 0.5 - L / 2, 0.5 - L / 2], 2), 2)
        run_case("Test8", rve, (Cut(geo7, geo4, 2),)) do
            v = group_volume(2)
            @test 0 < v < L^3
        end
    end

    @testset "Non-unit, anisotropic, shifted cell" begin
        # regression: setPeriodic used a unit translation regardless of the cell size
        rve = RVE([2.0, 1.0, 0.5], [1, 1, 1], [1.0, -3.0, 0.25], 0.15)
        incs = (Sphere([3.0, -2.5, 0.5], 0.2, 2),)  # crosses the xmax face
        run_case("TestNonUnit", rve, incs) do
            @test group_volume(2) ≈ sphere_vol(0.2) rtol = 1e-6
        end
    end

    @testset "Aqua" begin
        Aqua.test_all(RVEs)
    end
end
