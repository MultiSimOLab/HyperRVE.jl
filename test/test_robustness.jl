# Regression tests from the review of the v2 rewrite.

@testset "small inclusion with a coarse global size (no refine)" begin
    for (r, h) in ((0.03, nothing), (0.05, 0.2), (0.1, 0.25), (0.2, 0.5))
        rve = RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), r)))
        generate(rve; h=h) do res
            @test res.stats.nelements > 0
            @test mesh_volume() ≈ 1 rtol = 1e-2
            @test all(c -> c.ok, res.periodicity)
        end
    end
end

@testset "barely crossing a periodic face" begin
    for δ in (1e-5, 1e-4, 1e-3), (c, s) in (((0.4, 0.5, 0.5), -1), ((0.6, 0.5, 0.5), 1),
                                           ((0.5, 0.4, 0.5), -1), ((0.5, 0.5, 0.6), 1))
        d = findfirst(x -> x != 0.5, c)
        R = 0.4 + δ
        res = generate(RVE(unitcell(), Phase(:p, Sphere(c, R))); h=0.15)
        @test volume(res, :p) ≈ sphere_vol(R) rtol = 1e-6
        @test all(x -> x.ok, res.periodicity)
    end
end

@testset "bodies flush with or tangent to a periodic face" begin
    flush = [Box((0, 0, 0), (0.5, 0.5, 0.5)), Cylinder((0.5, 0.5, 0.5), (0.5, 0, 0), 0.2),
             Box((0.25, 0.25, 0.5), (0.5, 0.5, 0.5))]
    for g in flush
        res = generate(RVE(unitcell(), Phase(:p, g)); h=0.15)
        @test all(x -> x.ok, res.periodicity)
    end
    # a sphere tangent to a face from inside: either a valid periodic mesh or a PeriodicityError
    ok = try
        all(x -> x.ok, generate(RVE(unitcell(), Phase(:p, Sphere((0.75, 0.5, 0.5), 0.25))); h=0.15).periodicity)
    catch e
        e isa PeriodicityError
    end
    @test ok
end

@testset "meshed at small and large scales, no matrix, overlap with an image" begin
    for k in (1e-3, 1e3)
        cell = Cell((k, 2k, k); origin=(3k, -k, 0))
        res = generate(RVE(cell, Phase(:p, Sphere(cellpoint(cell, (0, 0, 0)), 0.3k))); h=0.2k)
        @test all(x -> x.ok, res.periodicity)
        @test volume_fraction(res, :p) ≈ sphere_vol(0.3) / 2 rtol = 1e-6
        @test res.stats.hmax < k
    end
    lattice = RVE(unitcell(), Phase(:solid, Sphere((0.5, 0.5, 0.5), 0.6)); matrix=nothing)
    generate(lattice; h=0.15) do res
        @test all(x -> x.ok, res.periodicity)
        @test mesh_volume() ≈ volume(res, :solid) rtol = 2e-2
    end
    # phase b overlaps only the periodic image of phase a
    @test_throws ErrorException generate(RVE(unitcell(), [Phase(:a, Sphere((0.95, 0.5, 0.5), 0.2)),
                                                          Phase(:b, Sphere((0.1, 0.5, 0.5), 0.1))];
                                             overlaps=:error); geometry_only=true)
end

@testset "verification detects non-periodic meshes" begin
    per = RVE(unitcell(), Phase(:p, Sphere((0, 0.5, 0.5), 0.25)))
    generate(RVE(unitcell(periodic=false), Phase(:p, Sphere((0, 0.5, 0.5), 0.25))); h=0.15) do res
        @test_throws PeriodicityError RVEs.verify_periodicity(per, RVEs.Frame(per))
    end
end

@testset "sessions: repeated keep_open, caller options, files" begin
    rve1 = RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.2)))
    rve2 = RVE(unitcell(), Phase(:q, Sphere((0.5, 0.5, 0.5), 0.3)))
    r1 = generate(rve1; h=0.3, keep_open=true)
    r2 = generate(rve2; h=0.1, keep_open=true)
    @test r1.model != r2.model
    f = joinpath(OUTDIR, "first.msh")
    write_mesh(r1, f)
    close_gmsh()
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f)
        @test length(gmsh.model.mesh.getNodes()[1]) == r1.stats.nnodes
        names = [gmsh.model.getPhysicalName(3, t) for (_, t) in gmsh.model.getPhysicalGroups(3)]
        @test "p" in names && !("q" in names)
    finally
        gmsh.finalize()
    end
    # caller-owned session: options restored
    gmsh.initialize(String[], false)
    try
        opts = ("General.Terminal", "Mesh.MeshSizeFromPoints", "Mesh.MshFileVersion", "Geometry.OCCBoundsUseStl")
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        before = [gmsh.option.getNumber(o) for o in opts]
        generate(rve1; h=0.3)
        @test [gmsh.option.getNumber(o) for o in opts] == before
    finally
        gmsh.finalize()
    end
    # written files are MSH 4.1 with periodic links, whatever the user configuration says
    f2 = joinpath(OUTDIR, "per41.msh")
    generate(rve1; h=0.3, output=f2)
    @test startswith(readlines(f2)[2], "4.1")
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f2)
        linked = count(gmsh.model.getEntities(2)) do (d, t)
            gmsh.model.mesh.getPeriodicNodes(d, t)[1] != t
        end
        @test linked > 0
    finally
        gmsh.finalize()
    end
    # geometry_only with output warns; do-block with keep_open keeps the session
    @test_logs (:warn, r"no mesh") generate(rve1; geometry_only=true, output=joinpath(OUTDIR, "geo.msh"))
    generate(rve1; geometry_only=true, keep_open=true) do r end
    @test gmsh.isInitialized() == 1
    close_gmsh()
end

@testset "known limitation: refinement close to a periodic face" begin
    # the size field is not periodic: a refined inclusion close to (but not crossing) a face
    # gets a fine mesh on one side only; the periodic copy of the coarse face degrades quality
    res = generate(RVE(unitcell(), Phase(:p, Sphere((0.78, 0.5, 0.5), 0.2); refine=InterfaceRefinement(h=0.02))); h=0.15)
    @test_broken res.stats.min_quality > 1e-2
end
