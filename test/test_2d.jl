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
    @test HyperRVE.characteristic_length(Ellipse((0, 0), (0.3, 0.1))) == 0.1
    @test HyperRVE.refpoint(Rectangle((0, 0), (2, 4))) == (1.0, 2.0)
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
        cellgroups(dim) = Set(filter(!startswith("interface"), group_names(res; dim=dim)))
        @test cellgroups(1) == Set(["xmin", "xmax", "ymin", "ymax", "boundary"])
        corners = ["corner_xmin_ymin", "corner_xmax_ymin", "corner_xmin_ymax", "corner_xmax_ymax"]
        @test cellgroups(0) == Set(vcat(["xmin", "xmax", "ymin", "ymax", "boundary", "corners"], corners))
        @test [group_tag(res, c; dim=0) for c in corners] == collect(31:34)
        @test [group_tag(res, f; dim=1) for f in ("xmin", "xmax", "ymin", "ymax")] == [1, 2, 3, 4]
        @test !any(startswith("edge_"), group_names(res))
        g(name, dim) = only(x for x in res.groups if x.name == name && x.dim == dim)
        @test length(g("corners", 0).entities) == 4
        s = HyperRVE.Frame(res.rve).scale
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
        @test mesh_area() ≈ 2.0 rtol = 1e-9          # straight outer boundary: exact
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

"Bounding box (x0, y0, x1, y1) and centre of mass of the geometry of a phase (unit cell: frame = identity)."
function phase_box(res, name)
    ents = res[name].entities
    b = [gmsh.model.occ.getBoundingBox(2, t) for t in ents]
    x0 = minimum(x[1] for x in b); y0 = minimum(x[2] for x in b)
    x1 = maximum(x[4] for x in b); y1 = maximum(x[5] for x in b)
    m = [gmsh.model.occ.getMass(2, t) for t in ents]
    c = reduce((a, b) -> a .+ b, (m[i] .* gmsh.model.occ.getCenterOfMass(2, ents[i]) for i in eachindex(ents))) ./ sum(m)
    return (x0, y0, x1, y1), c
end

@testset "positions and orientations (T1, T2)" begin
    sq = square()
    generate(RVE(sq, [Phase(:tall, Ellipse((0.3, 0.5), (0.1, 0.3))),             # rx < ry: tall
                      Phase(:wide, Ellipse((0.75, 0.2), (0.2, 0.05)))]); geometry_only=true) do res
        (b, _) = phase_box(res, :tall)
        @test b[3] - b[1] ≈ 0.2 atol = 1e-3
        @test b[4] - b[2] ≈ 0.6 atol = 1e-3
        (b, _) = phase_box(res, :wide)
        @test b[3] - b[1] ≈ 0.4 atol = 1e-3
        @test b[4] - b[2] ≈ 0.1 atol = 1e-3
    end
    rot = Rotate(Disk((0.7, 0.5), 0.1); angle=π / 2, center=(0.5, 0.5))           # -> (0.5, 0.7)
    cw = Rotate(Disk((0.7, 0.5), 0.1); angle=π / 2, center=(0.5, 0.5), axis=(0, 0, -1))  # -> (0.5, 0.3)
    tr = Translate(Disk((0.2, 0.2), 0.05), (0.3, 0.1))                            # -> (0.5, 0.3) ...
    sc = Scale(Disk((0.65, 0.35), 0.05), 2; center=(0.5, 0.5))                     # -> (0.8, 0.2), r 0.1
    generate(RVE(sq, [Phase(:rot, rot), Phase(:cw, cw), Phase(:tr, Translate(tr, (-0.3, 0.5))), Phase(:sc, sc)]);
             geometry_only=true) do res
        @test all(isapprox.(phase_box(res, :rot)[2][1:2], (0.5, 0.7); atol=1e-6))
        @test all(isapprox.(phase_box(res, :cw)[2][1:2], (0.5, 0.3); atol=1e-6))
        @test all(isapprox.(phase_box(res, :tr)[2][1:2], (0.2, 0.8); atol=1e-6))
        @test all(isapprox.(phase_box(res, :sc)[2][1:2], (0.8, 0.2); atol=5e-3))   # dilation: B-spline, OCC integrals ~1e-3
        @test volume(res, :sc) ≈ disk_area(0.1) rtol = 1e-2
    end
    @test_throws ArgumentError Rotate(Disk((0, 0), 1); angle=1.0, axis=(1, 0, 0))
    @test_throws ArgumentError Disk((0, 0), 1) ∪ Sphere((0, 0, 0), 1)
    @test_throws ArgumentError Disk((0, 0), 1) ∩ Sphere((0, 0, 0), 1)
    @test_throws ArgumentError setdiff(Disk((0, 0), 1), Sphere((0, 0, 0), 1))
end

@testset "side and corner positions (T3)" begin
    cell = Cell((2.0, 1.0); origin=(1.0, -1.0))
    generate(RVE(cell, Phase(:p, Disk(cellpoint(cell, (0.5, 0.5)), 0.2))); geometry_only=true) do res
        fr = HyperRVE.Frame(res.rve)
        phys(x) = fr.origin .+ fr.scale .* x
        g(name, dim) = only(x for x in res.groups if x.name == name && x.dim == dim)
        for (name, d, val) in (("xmin", 1, 1.0), ("xmax", 1, 3.0), ("ymin", 2, -1.0), ("ymax", 2, 0.0))
            for t in g(name, 1).entities
                b = gmsh.model.getBoundingBox(1, t)
                @test isapprox(phys((b[1], b[2]))[d], val; atol=1e-6)
                @test isapprox(phys((b[4], b[5]))[d], val; atol=1e-6)
            end
        end
        for (name, p) in (("corner_xmin_ymin", (1.0, -1.0)), ("corner_xmax_ymin", (3.0, -1.0)),
                          ("corner_xmin_ymax", (1.0, 0.0)), ("corner_xmax_ymax", (3.0, 0.0)))
            v = gmsh.model.getValue(0, only(g(name, 0).entities), Float64[])
            @test all(isapprox.(phys((v[1], v[2])), p; atol=1e-9))
        end
    end
end

@testset "reopened 2D file (T4)" begin
    f = joinpath(OUTDIR, "cell2d_full.msh")
    generate(RVE(Cell((2.0, 1.0); origin=(1.0, -1.0)), Phase(:p, Disk((1.0, -1.0), 0.3))); h=0.2, output=f)
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f)
        x = reshape(gmsh.model.mesh.getNodes()[2], 3, :)
        @test all(iszero, x[3, :])
        @test extrema(x[1, :]) == (1.0, 3.0) || all(extrema(x[1, :]) .≈ (1.0, 3.0))
        @test all(extrema(x[2, :]) .≈ (-1.0, 0.0))
        @test isempty(gmsh.model.mesh.getElements(3)[1])
        names(dim) = Set(gmsh.model.getPhysicalName(d, t) for (d, t) in gmsh.model.getPhysicalGroups(dim))
        @test names(2) == Set(["matrix", "p"])
        @test Set(["xmin", "xmax", "ymin", "ymax", "boundary"]) ⊆ names(1)
        @test "corner_xmax_ymax" in names(0)
    finally
        gmsh.finalize()
    end
end

@testset "more 2D cases (T5-T9)" begin
    sq = square()
    # no matrix: porous / lattice-like
    res = generate(RVE(sq, Phase(:solid, Disk((0.5, 0.5), 0.6)); matrix=nothing); h=0.1)
    @test phase_names(res) == ["solid"]
    @test all(c -> c.ok, res.periodicity)
    # partial and no periodicity
    part = generate(RVE(Cell((1.0, 1.0); periodic=(true, false)), Phase(:p, Disk((0.5, 0.0), 0.2))); h=0.1)
    @test [c.direction for c in part.periodicity] == [1]
    @test volume(part, :p) ≈ disk_area(0.2) / 2 rtol = 1e-9
    nop = generate(RVE(Cell((1.0, 1.0); periodic=false), Phase(:p, Disk((0.5, 0.5), 0.2))); h=0.1)
    @test isempty(nop.periodicity) && nop.stats.nelements > 0
    # tangent disks, laminate, rectangle = cell
    for g in (Disk((0.5, 0.5), 0.5), Rectangle((0.0, 0.0), (1.0, 0.4)), Disk((0.75, 0.5), 0.25))
        r = generate(RVE(sq, Phase(:p, g)); h=0.1)
        @test all(c -> c.ok && c.npairs > 0, r.periodicity)
    end
    full = @test_logs (:warn, r"matrix") generate(RVE(sq, Phase(:p, Rectangle((0, 0), (1, 1)))); h=0.2)
    @test volume_fraction(full, :p) ≈ 1
    # refinement and curvature act in 2D
    base = generate(RVE(sq, Phase(:p, Disk((0.5, 0.5), 0.2))); h=0.1).stats.nelements
    @test generate(RVE(sq, Phase(:p, Disk((0.5, 0.5), 0.2); refine=InterfaceRefinement(level=3))); h=0.1).stats.nelements > 2base
    @test generate(RVE(sq, Phase(:p, Disk((0.5, 0.5), 0.02))); h=0.1).stats.hmin < 0.03
    # booleans and overlaps in 2D
    @test volume(geo(RVE(sq, Phase(:p, Fuse(Disk((0.4, 0.5), 0.1), Disk((0.6, 0.5), 0.1))))), :p) ≈ 2 * disk_area(0.1) rtol = 1e-9
    @test volume(geo(RVE(sq, Phase(:p, Intersect(Rectangle((0, 0), (0.5, 1)), Disk((0.5, 0.5), 0.2))))), :p) ≈ disk_area(0.2) / 2 rtol = 1e-9
    @test_throws ErrorException geo(RVE(sq, [Phase(:a, Disk((0.5, 0.5), 0.2)), Phase(:b, Disk((0.6, 0.5), 0.2))]; overlaps=:error))
end
