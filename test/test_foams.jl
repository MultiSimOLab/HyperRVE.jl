# Void phases, interface groups as a union, closed- and open-cell foams.

@testset "void phases" begin
    c = unitcell()
    res = generate(RVE(c, Phase(:pore, Sphere((0.5, 0.5, 0.5), 0.3); void=true)); h=0.1)
    @test [p.name for p in res.phases] == ["matrix"]
    @test res.phases[1].fraction ≈ 1 - 4π * 0.027 / 3 rtol = 1e-6
    @test interface_names(res) == ["interface_matrix_void"]
    # a void phase before other phases; across a periodic face
    res = generate(RVE(c, [Phase(:pore, Sphere((0.05, 0.5, 0.5), 0.15); void=true), Phase(:inc, Sphere((0.6, 0.6, 0.6), 0.2))]); h=0.1)
    @test [p.name for p in res.phases] == ["matrix", "inc"]
    @test sort(interface_names(res)) == ["interface_matrix_inc", "interface_matrix_void"]
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume(res, :inc) ≈ 4π * 0.008 / 3 rtol = 1e-6
    # painted over an earlier phase
    g = geo(RVE(c, [Phase(:inc, Sphere((0.5, 0.5, 0.5), 0.3)), Phase(:hole, Sphere((0.5, 0.5, 0.5), 0.1); void=true)]))
    @test volume(g, :inc) ≈ 4π * (0.027 - 0.001) / 3 rtol = 1e-6
    @test_throws ArgumentError Phase(:p, Sphere((0.5, 0.5, 0.5), 0.1); void=true, h=0.01)
    @test_throws ArgumentError Phase(:p, Sphere((0.5, 0.5, 0.5), 0.1); void=true, refine=0.01)
    # 2D
    c2 = Cell((1.0, 1.0))
    res = generate(RVE(c2, Phase(:hole, Disk((0.5, 0.5), 0.2); void=true)); h=0.05)
    @test res.phases[1].fraction ≈ 1 - π * 0.04 rtol = 1e-4
end

@testset "interfaces = :union" begin
    t = voronoi(unitcell(); n=6, rng=1)
    res = generate(RVE(unitcell(), grain_phases(t); matrix=nothing); h=0.15, interfaces=:union)
    names = unique(g.name for g in res.groups)
    @test "interfaces" in names && !any(n -> startswith(n, "interface_"), names)
    @test length(res.interfaces) > 6 && all(i -> i.tag == 0, res.interfaces)
    @test_throws ArgumentError generate(RVE(unitcell(), Phase(:s, Sphere((0.5, 0.5, 0.5), 0.2))); geometry_only=true, interfaces=:pairs)
end

@testset "foams" begin
    c = unitcell()
    t = voronoi(c; n=8, rng=1)
    # closed cells: exact density, empty pores, periodic
    w = wall_thickness(t, 0.15)
    p = pores(t; wall=w)
    @test 1 - sum(abs(HyperRVE._polyhedron_volume(x.vertices, x.faces)) for x in p) ≈ 0.15 atol = 1e-9
    res = generate(RVE(c, Phase(:pores, p; void=true); matrix=:walls); h=w)
    @test res.phases[1].fraction ≈ 0.15 atol = 1e-8
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test isempty(pores(t; wall=0.0)) == false && length(pores(t; wall=0.0)) == 8
    @test_logs (:warn, r"vanish") match_mode = :any pores(t; wall=0.9)
    @test_throws ArgumentError wall_thickness(t, 1.2)
    # 2D closed foam
    t2 = voronoi(Cell((1.0, 1.0)); n=10, rng=2)
    w2 = wall_thickness(t2, 0.2)
    res = generate(RVE(t2.cell, Phase(:pores, pores(t2; wall=w2); void=true); matrix=:walls); h=w2 / 2)
    @test res.phases[1].fraction ≈ 0.2 atol = 1e-8
    # open cells: the edges as struts
    nodes, edges = tessellation_edges(t)
    @test length(nodes) - length(edges) + sum(length, t.loops) ÷ 2 - 8 == 0       # Euler on the 3-torus
    t4 = voronoi(c; n=4, rng=2)
    its = foam_struts(t4; radius=0.03)
    @test count(x -> x isa Cylinder, its) <= length(tessellation_edges(t4)[2])
    g = geo(RVE(c, Phase(:foam, its; fuse=true); matrix=nothing))
    @test 0.02 < volume_fraction(g, :foam) < relative_density(c, tessellation_edges(t4)...; radius=0.03)
    @test_throws ArgumentError foam_struts(voronoi(Cell((1.0, 1.0, 1.0); periodic=false); n=4, rng=1); radius=0.03)
end

@testset "foam_struts keeps the radius (closure regression)" begin
    t = voronoi(unitcell(); n=8, lloyd=2, rng=1)
    its = foam_struts(t; radius=0.025)
    @test all(x -> x.radius == 0.025, filter(x -> x isa Cylinder, its))
end

@testset "foams and voids: review fixes" begin
    # a block cut from a periodic tessellation: pores clipped and wrapped, exact density
    cn = Cell((1.0, 1.0, 1.0); periodic=false)
    t = voronoi(cn; n=6, rng=1, periodic=true)
    w = wall_thickness(t, 0.2)
    g = geo(RVE(cn, Phase(:pores, pores(t; wall=w); void=true); matrix=:walls))
    @test g.phases[1].fraction ≈ 0.2 atol = 1e-8
    # everything void and no matrix: nothing to mesh
    @test_throws ArgumentError RVE(unitcell(), Phase(:pore, Sphere((0.5, 0.5, 0.5), 0.3); void=true); matrix=nothing)
    # a void slab covering a pair of periodic faces is not a periodicity error
    res = generate(RVE(unitcell(), Phase(:slab, Box((-0.1, -1.0, -1.0), (0.2, 3.0, 3.0)); void=true)); h=0.1)
    @test res.phases[1].fraction ≈ 0.8 rtol = 1e-8
    # nearly vanished pores are skipped, not an error
    tc = voronoi(unitcell(); seeds=[(x, y, z) for x in (0.25, 0.75) for y in (0.25, 0.75) for z in (0.25, 0.75)])
    @test length(pores(tc; wall=0.5 - 1e-9)) == 0
    @test all(p -> all(isapprox.(abs(HyperRVE._polyhedron_volume(p.vertices, p.faces)), 0.4^3)), pores(tc; wall=0.1))
    # interfaces keyword validated before meshing; interface_tag with :union
    @test_throws ArgumentError generate(RVE(unitcell(), Phase(:s, Sphere((0.5, 0.5, 0.5), 0.2))); geometry_only=true, interfaces=1)
    r = generate(RVE(unitcell(), Phase(:s, Sphere((0.5, 0.5, 0.5), 0.2))); geometry_only=true, interfaces=:union)
    @test_throws ErrorException interface_tag(r, :matrix, :s)
    @test_logs (:warn, r"min_length") match_mode = :any foam_struts(voronoi(unitcell(); n=2, rng=1); radius=0.02, min_length=0.45)
end
