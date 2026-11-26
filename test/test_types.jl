using RVEs: refpoint, characteristic_length, isunbounded, invariant_axis, _shift_range

@testset "constructors and conversion" begin
    c = Cell([1, 2, 3])
    @test c isa Cell{3}
    @test c.size === (1.0, 2.0, 3.0)
    @test c.origin === (0.0, 0.0, 0.0)
    @test c.periodic === (true, true, true)
    @test Cell(size=(1, 1, 1), periodic=(true, false, true)).periodic === (true, false, true)
    @test Cell((1.0, 1.0)) isa Cell{2}
    @test cellpoint(Cell((2.0, 1.0, 1.0); origin=(1, 1, 1)), (0.5, 0.5, 0)) == (2.0, 1.5, 1.0)
    @test volume(Cell((2.0, 1.0, 0.5))) == 1.0
    @test Sphere([0, 0, 0], 1) == Sphere(center=(0, 0, 0), radius=1.0)
    @test Box(center=(0.5, 0.5, 0.5), size=(1, 1, 1)).corner == (0.0, 0.0, 0.0)
    @test Box(corner=(0, 0, 0), size=(1, 2, 3)) == Box((0, 0, 0), (1, 2, 3))
    cyl = Cylinder(center=(0, 0, 0), direction=(0, 0, 2), length=4, radius=1)
    @test cyl.base == (0.0, 0.0, -2.0) && cyl.axis == (0.0, 0.0, 4.0)
    @test Fibre((0, 0, 0), (0, 3, 0), 0.1).direction == (0.0, 1.0, 0.0)
    s = Sphere((0, 0, 0), 1)
    @test Fuse(s, s) isa Fuse{3}
    @test Fuse([s, s]).parts == Fuse(s, s).parts
    @test s ∪ s isa Fuse{3}
    @test union(s, s, s) isa Fuse{3} && length(union(s, s, s).parts) == 3
    @test length(intersect(s, s, s).parts) == 3
    @test s ∩ s isa Intersect{3}
    @test setdiff(s, s) isa Cut{3}
    @test Scale(s, 2).factors == (2.0, 2.0, 2.0)
    p = Phase("p", (Sphere((0, 0, 0), 0.1) for _ in 1:3))
    @test p isa Phase{3} && length(p.items) == 3
    @test Phase(:q, s).name == "q"
    @test MeshOptions(h=0.1).h == 0.1
end

@testset "validation" begin
    @test_throws ArgumentError Sphere((0, 0, 0), -1)
    @test_throws ArgumentError Sphere((0, 0), 1)
    @test_throws ArgumentError Sphere((0, 0, NaN), 1)
    @test_throws ArgumentError Box((0, 0, 0), (1, 0, 1))
    @test_throws ArgumentError Box(size=(1, 1, 1))
    @test_throws ArgumentError Box(corner=(0, 0, 0), center=(0, 0, 0), size=(1, 1, 1))
    @test_throws ArgumentError Cylinder((0, 0, 0), (0, 0, 0), 1)
    @test_throws ArgumentError Fibre((0, 0, 0), (0, 0, 0), 1)
    @test_throws ArgumentError Ellipsoid((0, 0, 0), (1, -1, 1))
    @test_throws ArgumentError Cell((1.0, 0.0, 1.0))
    @test_throws ArgumentError Cell((1.0,))
    @test_throws ArgumentError Cell((1.0, 1.0, 1.0); periodic=(true, false))
    @test_throws ArgumentError Fuse(Geometry{3}[])
    @test_throws ArgumentError Scale(Sphere((0, 0, 0), 1), 0)
    @test_throws ArgumentError Phase("a b", Sphere((0, 0, 0), 1))
    @test_throws ArgumentError Phase("", Sphere((0, 0, 0), 1))
    @test_throws ArgumentError InterfaceRefinement(level=0)
    @test_throws ArgumentError InterfaceRefinement(distmin=0.5, distmax=0.1)
    s = Sphere((0.5, 0.5, 0.5), 0.1)
    c = unitcell()
    @test_throws ArgumentError RVE(c, [Phase(:a, s), Phase(:a, s)])
    @test_throws ArgumentError RVE(c, Phase(:matrix, s))
    @test_throws ArgumentError RVE(c, Phase(:xmin, s))
    @test_throws ArgumentError RVE(c, Phase(:corner_1, s))
    @test_throws ArgumentError RVE(c, Phase(:a, s); tiling=0)
    @test_throws ArgumentError RVE(c, Phase(:a, s); overlaps=:maybe)
    # unbounded, not axis-aligned
    @test_throws ArgumentError RVE(c, Phase(:f, Fibre((0.5, 0.5, 0.5), (1, 1, 0), 0.1)))
    cross = Fuse(Fibre((0.5, 0.5, 0.5), (1, 0, 0), 0.1), Fibre((0.5, 0.5, 0.5), (0, 1, 0), 0.1))
    @test_throws ArgumentError RVE(c, Phase(:f, cross))
    # not implemented yet in this version
    @test_throws ArgumentError RVE(Cell((1.0, 1.0)), Phase(:a, s))
    @test_throws ArgumentError RVE(c, Phase(:a, s); tiling=2)
    @test_throws ArgumentError generate(RVE(c, Phase(:a, s)); mesh=MeshOptions(order=2))
    @test_throws ArgumentError generate(RVE(c, Phase(:a, s)); mesh=MeshOptions(optimize=false))
    @test_throws ArgumentError generate(RVE(c, Phase(:a, s)); output=joinpath(OUTDIR, "nodir", "x.msh"))
    @test_throws ArgumentError generate(RVE(c, Phase(:a, s)); h=0.1, mesh=MeshOptions(h=0.2))
    @test_logs (:warn, r"Fibre") RVE(c, Phase(:a, Cylinder((-0.1, 0.5, 0.5), (2, 0, 0), 0.1)))
end

@testset "traits" begin
    s = Sphere((1, 2, 3), 0.5)
    @test refpoint(s) == (1.0, 2.0, 3.0)
    @test refpoint(Box((0, 0, 0), (2, 2, 2))) == (1.0, 1.0, 1.0)
    @test refpoint(Cylinder((0, 0, 0), (0, 0, 2), 1)) == (0.0, 0.0, 1.0)
    @test refpoint(Translate(s, (1, 0, 0))) == (2.0, 2.0, 3.0)
    @test all(refpoint(Rotate(s; angle=π / 2, axis=(0, 0, 1), center=(0, 0, 0))) .≈ (-2.0, 1.0, 3.0))
    @test refpoint(Rotate(s; angle=1.0)) == s.center      # own centre by default
    @test refpoint(Scale(s, 2; center=(0, 0, 0))) == (2.0, 4.0, 6.0)
    @test characteristic_length(Ellipsoid((0, 0, 0), (3, 1, 2))) == 1.0
    @test characteristic_length(Cut(Box((0, 0, 0), (1, 1, 1)), Sphere((0, 0, 0), 0.2))) == 0.2
    @test characteristic_length(Scale(s, (2, 3, 4))) == 1.0
    f = Fibre((0, 0, 0), (0, 1, 0), 0.1)
    @test isunbounded(f) && invariant_axis(f) == 2
    @test !isunbounded(s) && invariant_axis(s) == 0
    @test isunbounded(Fuse(f, s)) && invariant_axis(Fuse(f, s)) == 0
    @test !isunbounded(Intersect(f, s))
    @test invariant_axis(Cut(f, Fibre((0, 0, 0), (0, 1, 0), 0.05))) == 2
    @test invariant_axis(Cut(f, s)) == 0
    @test invariant_axis(Rotate(f; angle=0.3, axis=(0, 1, 0))) == 2
    @test invariant_axis(Rotate(f; angle=0.3, axis=(0, 0, 1))) == 0
    @test invariant_axis(Translate(f, (1, 1, 1))) == 2
end

@testset "shift ranges" begin
    t = 1e-6
    @test _shift_range(0.3, 0.7, 1.0, 1, true, false, t) == [0]            # interior
    @test _shift_range(-0.2, 0.2, 1.0, 1, true, false, t) == [0, 1]        # crosses min
    @test _shift_range(0.8, 1.2, 1.0, 1, true, false, t) == [-1, 0]        # crosses max
    @test _shift_range(-0.1, 1.1, 1.0, 1, true, false, t) == [-1, 0, 1]    # larger than cell
    @test _shift_range(-0.8, 1.8, 1.0, 1, true, false, t) == [-1, 0, 1]
    @test _shift_range(0.6, 1.0, 1.0, 1, true, false, t) == [-1, 0]        # flush: image imprints xmin
    @test _shift_range(0.0, 0.4, 1.0, 1, true, false, t) == [0, 1]
    @test _shift_range(-0.2, 0.2, 1.0, 1, false, false, t) == [0]          # non-periodic
    @test _shift_range(1.2, 1.4, 1.0, 1, false, false, t) == Int[]         # outside
    @test _shift_range(-0.2, 0.2, 1.0, 1, true, true, t) == 0:0            # invariant
    @test _shift_range(0.8, 1.2, 1.0, 3, true, false, t) == [-1, 0, 1, 2]  # tiling
    @test _shift_range(0.8, 1.2, 1.0, 3, false, false, t) == [0, 1, 2]
    @test _shift_range(2.8, 3.2, 1.0, 3, false, false, t) == [0]
end

@testset "show" begin
    s = Sphere((0.5, 0.5, 0.5), 0.1)
    r = RVE(unitcell(), Phase(:p, s; refine=InterfaceRefinement()))
    @test occursin("p", sprint(show, r))
    @test occursin("periodic xyz", sprint(show, unitcell()))
    @test occursin("refined", sprint(show, r.phases[1]))
end
