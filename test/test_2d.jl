# 2D cells: disks, ellipses, rectangles; groups on curves and points; triangles.

disk_area(r) = π * r^2
square() = Cell((1.0, 1.0))

@testset "constructors and validation" begin
    @test Disk((0, 0), 1) == Disk(center=(0, 0), radius=1)
    @test Rectangle(center=(0.5, 0.5), size=(1, 1)).corner == (0.0, 0.0)
    @test_throws ArgumentError Disk((0, 0, 0), 1)
    @test_throws ArgumentError Ellipse((0, 0), (1, 0))
    @test_throws ArgumentError Rectangle(size=(1, 1))
    @test_throws ArgumentError Phase(:p, [Disk((0, 0), 1), Sphere((0, 0, 0), 1)])
    @test_throws ArgumentError RVE(square(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.1)))
    @test_throws ArgumentError RVE(Cell((1.0, 1.0, 1.0)), Phase(:p, Disk((0.5, 0.5), 0.1)))
    @test RVEs.characteristic_length(Ellipse((0, 0), (0.3, 0.1))) == 0.1
    @test RVEs.refpoint(Rectangle((0, 0), (2, 4))) == (1.0, 2.0)
end

@testset "areas (B1)" begin
    r = 0.2
    for c in ((0.5, 0.5), (0.0, 0.5), (0.0, 0.0), (1.0, 1.0))
        @test volume(geo(RVE(square(), Phase(:p, Disk(c, r)))), :p) ≈ disk_area(r) rtol = 1e-9
    end
    hex = Cell((1.0, sqrt(3.0)))
    res = geo(RVE(hex, Phase(:fibre, [Disk(cellpoint(hex, (0.5, 0.5)), 0.3), Disk(cellpoint(hex, (0, 0)), 0.3)]);
                  matrix=:resin))
    @test volume_fraction(res, :fibre) ≈ 2π * 0.09 / sqrt(3) rtol = 1e-9
    @test phase_names(res) == ["resin", "fibre"]
    # ellipses with either radius larger, rotated, crossing a corner
    for radii in ((0.3, 0.1), (0.1, 0.3))
        e = Rotate(Ellipse((0.0, 0.0), radii); angle=0.4)
        @test volume(geo(RVE(square(), Phase(:e, e))), :e) ≈ π * 0.03 rtol = 1e-6
    end
    @test volume(geo(RVE(square(), Phase(:p, Cut(Rectangle(center=(0.5, 0.5), size=(0.6, 0.6)),
                                                    Disk((0.5, 0.5), 0.2))))), :p) ≈ 0.36 - disk_area(0.2) rtol = 1e-9
    # non-uniform scaling gives a B-spline curve: OCC areas accurate to ~1e-2
    @test volume(geo(RVE(square(), Phase(:p, Scale(Disk((0.5, 0.5), 0.1), (2, 1))))), :p) ≈ 2 * disk_area(0.1) rtol = 1e-2
    @test volume(geo(RVE(square(), Phase(:p, Translate(Disk((0.9, 0.5), 0.2), (1, 0))))), :p) ≈ disk_area(0.2) rtol = 1e-9
    # flush with a side
    @test volume(geo(RVE(square(), Phase(:p, Rectangle((0, 0), (0.5, 0.5))))), :p) ≈ 0.25 rtol = 1e-12
end

@testset "groups (B2)" begin
    cell = Cell((2.0, 1.0); origin=(1.0, -1.0))
    generate(RVE(cell, Phase(:p, Disk(cellpoint(cell, (0, 0)), 0.2))); geometry_only=true) do res
        @test sort(group_names(res; dim=2)) == ["matrix", "p"]
        @test Set(group_names(res; dim=1)) == Set(["xmin", "xmax", "ymin", "ymax", "boundary"])
        corners = ["corner_xmin_ymin", "corner_xmax_ymin", "corner_xmin_ymax", "corner_xmax_ymax"]
        @test Set(group_names(res; dim=0)) == Set(vcat(["xmin", "xmax", "ymin", "ymax", "boundary", "corners"], corners))
        @test [group_tag(res, c; dim=0) for c in corners] == collect(31:34)
        @test [group_tag(res, f; dim=1) for f in ("xmin", "xmax", "ymin", "ymax")] == [1, 2, 3, 4]
        @test !any(startswith("edge_"), group_names(res))
        g(name, dim) = only(x for x in res.groups if x.name == name && x.dim == dim)
        @test length(g("corners", 0).entities) == 4
        s = RVEs.Frame(res.rve).scale
        len(name) = s * sum(gmsh.model.occ.getMass(1, t) for t in g(name, 1).entities)
        @test len("xmin") ≈ 1.0 rtol = 1e-9
        @test len("ymax") ≈ 2.0 rtol = 1e-9
    end
end

@testset "periodic triangle meshes (B3)" begin
    cell = Cell((2.0, 1.0); origin=(1.0, -1.0))
    rve = RVE(cell, [Phase(:p, Disk(cellpoint(cell, (0, 0)), 0.2); refine=InterfaceRefinement(level=2)),
                     Phase(:q, Rotate(Ellipse(cellpoint(cell, (0.5, 1.0)), (0.3, 0.1)); angle=0.3))])
    generate(rve; h=0.1) do res
        @test length(res.periodicity) == 2
        @test all(c -> c.ok && c.npairs > 0, res.periodicity)
        @test brute_periodic_match_2d("xmin", "xmax", 1, 2.0)
        @test brute_periodic_match_2d("ymin", "ymax", 2, 1.0)
        @test mesh_area() ≈ 2.0 rtol = 1e-2
        @test res.stats.nelements > 0 && res.stats.min_quality > 0
        types, _, _ = gmsh.model.mesh.getElements(2)
        @test types == [2]                               # 3-node triangles only
        @test isempty(gmsh.model.mesh.getElements(3)[1])
        z = reshape(gmsh.model.mesh.getNodes()[2], 3, :)[3, :]
        @test all(iszero, z)
    end
    @test occursin("triangles", sprint(show, MIME("text/plain"), generate(rve; h=0.2)))
end

@testset "scale invariance and supercells (B4)" begin
    counts = Int[]
    for k in (1e-6, 1.0, 1e3)
        cell = Cell((k, k); origin=(2k, -k))
        res = generate(RVE(cell, Phase(:p, Disk(cellpoint(cell, (1, 0)), 0.3k))); h=0.1k)
        @test volume_fraction(res, :p) ≈ disk_area(0.3) rtol = 1e-9
        @test all(c -> c.ok && c.npairs > 0, res.periodicity)
        push!(counts, res.stats.nelements)
    end
    # small 2D meshes: the frame round-off changes a few triangles
    @test maximum(counts) - minimum(counts) <= 0.05 * minimum(counts)
    unit = RVE(square(), Phase(:p, [Disk((0.5, 0.5), 0.2), Disk((0, 0), 0.1)]))
    sup = generate(supercell(unit, (3, 2)); h=0.15)
    @test volume(sup) ≈ 6
    @test volume(sup, :p) ≈ 6 * (disk_area(0.2) + disk_area(0.1)) rtol = 1e-9
    @test all(c -> c.ok && c.npairs > 0, sup.periodicity)
    # file written in physical coordinates
    k = 1e-3
    f = joinpath(OUTDIR, "cell2d.msh")
    generate(RVE(Cell((k, 2k); origin=(5k, 0)), Phase(:p, Disk((5k, 0), 0.3k))); h=0.2k, output=f)
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f)
        x = reshape(gmsh.model.mesh.getNodes()[2], 3, :)
        @test minimum(x[1, :]) ≈ 5k && maximum(x[1, :]) ≈ 6k
        @test maximum(x[2, :]) ≈ 2k
    finally
        gmsh.finalize()
    end
end
