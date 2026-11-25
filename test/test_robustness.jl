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
        @test length(res[:p].entities) == 2     # the image was created
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    end
end

@testset "bodies flush with or tangent to a periodic face" begin
    flush = [Box((0, 0, 0), (0.5, 0.5, 0.5)), Cylinder((0.5, 0.5, 0.5), (0.5, 0, 0), 0.2),
             Box((0.25, 0.25, 0.5), (0.5, 0.5, 0.5))]
    vols = (0.125, π * 0.2^2 * 0.5, 0.125)
    for (g, v) in zip(flush, vols)
        res = generate(RVE(unitcell(), Phase(:p, g)); h=0.15)
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
        @test volume(res, :p) ≈ v rtol = 1e-6
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
    counts = Int[]
    for k in (1e-6, 1e-3, 0.1, 1.0, 1e3, 1e6)
        cell = Cell((k, 2k, k); origin=(3k, -k, 0.5k))
        res = generate(RVE(cell, Phase(:p, Sphere(cellpoint(cell, (0, 0, 0)), 0.3k))); h=0.2k)
        @test length(res.periodicity) == 3
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
        @test "xmin" in group_names(res; dim=2) && "corners" in group_names(res; dim=0)
        @test volume_fraction(res, :p) ≈ sphere_vol(0.3) / 2 rtol = 1e-6
        @test volume(res) ≈ 2k^3 rtol = 1e-12
        @test res.stats.hmax < k && res.stats.hmin > 0.001k
        push!(counts, res.stats.nelements)
    end
    # the same computation runs at every scale (normalised frame; the frame origin and scale
    # are not exactly representable, so allow tiny differences)
    @test maximum(counts) - minimum(counts) <= 0.01 * minimum(counts)
    # written files carry physical coordinates and periodic links
    k = 1e-6
    f = joinpath(OUTDIR, "micro.msh")
    generate(RVE(Cell((k, k, k); origin=(5k, 0, 0)), Phase(:p, Sphere((5k, 0, 0), 0.3k))); h=0.2k, output=f)
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f)
        x = reshape(gmsh.model.mesh.getNodes()[2], 3, :)
        @test minimum(x[1, :]) ≈ 5k rtol = 1e-9
        @test maximum(x[1, :]) ≈ 6k rtol = 1e-9
        linked = count(e -> gmsh.model.mesh.getPeriodicNodes(e...)[1] != e[2], gmsh.model.getEntities(2))
        @test linked > 0
    finally
        gmsh.finalize()
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
    # caller-owned session: options and current model restored; write_mesh still writes 4.1
    gmsh.initialize(String[], false)
    try
        opts = RVEs._SESSION_OPTIONS
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.model.add("mine")
        before = [gmsh.option.getNumber(o) for o in opts]
        r = generate(rve1; h=0.3)
        @test [gmsh.option.getNumber(o) for o in opts] == before
        @test gmsh.model.getCurrent() == "mine"
        @test_throws ErrorException generate(RVE(unitcell(), Phase(:p, Cut(Sphere((0.5, 0.5, 0.5), 0.1),
                                                                            Sphere((0.5, 0.5, 0.5), 0.2)))))
        @test gmsh.model.getCurrent() == "mine"
        f3 = joinpath(OUTDIR, "caller.msh")
        write_mesh(r, f3)
        @test startswith(readlines(f3)[2], "4.1")
        @test gmsh.option.getNumber("Mesh.MshFileVersion") == 2.2
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

@testset "B-spline bodies crossing periodic faces" begin
    bodies = [Ellipsoid((0.05, 0.5, 0.5), (0.2, 0.1, 0.15)),
              Rotate(Ellipsoid((0.5, 0.95, 0.5), (0.25, 0.1, 0.1)); axis=(1, 1, 1), angle=0.6),
              Scale(Sphere((0.5, 0.5, 0.02), 0.1), (2, 1, 1))]
    for g in bodies, m in (:matrix, nothing)
        res = generate(RVE(Cell((2.0, 1.0, 1.0); origin=(-1.0, 0.0, 0.0)), Phase(:p, Translate(g, (-1.0, 0.0, 0.0)));
                           matrix=m); h=0.12)
        @test Set(["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]) ⊆ Set(group_names(res; dim=2)) || m === nothing
        @test all(x -> x.ok, res.periodicity)
        crossing = findall(x -> x.npairs > 0, res.periodicity)
        @test !isempty(crossing)
    end
end

@testset "frame mapping of transforms in a shifted, scaled cell" begin
    phases(c, L) = [Phase(:e, Rotate(Rotate(Ellipsoid(c .+ (0.0, 0.0, 0.0), (0.25L, 0.1L, 0.1L)); axis=(0, 1, 0), angle=π / 4);
                                    axis=(0, 0, 1), angle=π / 4)),
                    Phase(:s, Scale(Sphere(c .+ (0.3L, 0.3L, 0.3L), 0.05L), (2, 1, 1))),
                    Phase(:t, Translate(Sphere(c, 0.08L), (0.4L, -0.4L, 0.1L))),
                    Phase(:b, Cut(Box(center=c .- (0.3L, 0.3L, 0.3L), size=(0.2L, 0.2L, 0.2L)),
                                  Fibre(c .- (0.3L, 0.3L, 0.3L), (1, 1, 1), 0.05L))),
                    Phase(:r, Rotate(Cylinder(center=c .+ (-0.3L, 0.3L, 0.0), direction=(0, 0, 1), length=0.3L, radius=0.05L);
                                     angle=0.7, axis=(1, 0, 0), center=c))]
    ref = geo(RVE(unitcell(), phases((0.5, 0.5, 0.5), 1.0)))
    for (L, o) in ((1e-3, (2.0, -1.0, 0.5)), (5e2, (1e3, 0.0, -7e2)))
        cell = Cell((L, L, L); origin=o)
        res = geo(RVE(cell, phases(cellpoint(cell, (0.5, 0.5, 0.5)), L)))
        for n in phase_names(ref)
            @test volume_fraction(res, n) ≈ volume_fraction(ref, n) rtol = 1e-8
        end
    end
end

@testset "msh output variants" begin
    k = 2e-6
    rve = RVE(Cell((k, k, k); origin=(3k, 0, 0)), Phase(:p, Sphere((3k, 0, 0), 0.3k)))
    f4 = joinpath(OUTDIR, "v.msh4")
    generate(rve; h=0.25k, output=f4) do r
        write_mesh(r, joinpath(OUTDIR, "w.msh"))
        write_mesh(r, joinpath(OUTDIR, "w.msh"))
        @test count(==(joinpath(OUTDIR, "w.msh")), r.files) == 1
    end
    for f in (f4, joinpath(OUTDIR, "w.msh"))
        gmsh.initialize(String[], false)
        try
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.open(f)
            x = reshape(gmsh.model.mesh.getNodes()[2], 3, :)
            @test maximum(x[1, :]) ≈ 4k rtol = 1e-9
            @test maximum(x[2, :]) ≈ k rtol = 1e-9
        finally
            gmsh.finalize()
        end
    end
    # caller session with binary output configured: write_mesh still writes ASCII 4.1
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Mesh.Binary", 1)
        r = generate(rve; h=0.25k)
        fb = joinpath(OUTDIR, "b.msh")
        write_mesh(r, fb)
        @test readlines(fb)[2] == "4.1 0 8"
        @test gmsh.option.getNumber("Mesh.Binary") == 1
    finally
        gmsh.finalize()
    end
end

@testset "curvature-based sizes only act on small radii" begin
    big = RVE(unitcell(), Phase(:p, [Sphere((0.5, 0.5, 0.5), 0.3)]))
    n0 = generate(big; h=0.1, mesh=MeshOptions(curvature=0)).stats.nelements
    n6 = generate(big; h=0.1).stats.nelements
    @test abs(n6 - n0) <= 0.01 * n0
    small = RVE(unitcell(), Phase(:p, [Sphere((0.5, 0.5, 0.5), 0.04)]))
    @test generate(small; h=0.1).stats.hmin < 0.05
end

@testset "known limitation: refinement close to a periodic face" begin
    # the size field is not periodic: a refined inclusion close to (but not crossing) a face
    # gets a fine mesh on one side only; the periodic copy of the coarse face degrades quality
    res = generate(RVE(unitcell(), Phase(:p, Sphere((0.78, 0.5, 0.5), 0.2); refine=InterfaceRefinement(h=0.02))); h=0.15)
    @test_broken res.stats.min_quality > 1e-2
end
