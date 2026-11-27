using HyperRVE, Test, Random
const HR = HyperRVE

lognormal_r(rng) = (while true
    d = 0.3 * exp(0.33 * randn(rng) - 0.05); 0.2 <= d <= 0.5 && return d / 2
end)

@testset "ball ∩ convex body volumes" begin
    box(lo, hi) = HR._Poly3([[(lo[1], lo[2], lo[3]), (lo[1], lo[2], hi[3]), (lo[1], hi[2], hi[3]), (lo[1], hi[2], lo[3])],
                             [(hi[1], lo[2], lo[3]), (hi[1], hi[2], lo[3]), (hi[1], hi[2], hi[3]), (hi[1], lo[2], hi[3])],
                             [(lo[1], lo[2], lo[3]), (hi[1], lo[2], lo[3]), (hi[1], lo[2], hi[3]), (lo[1], lo[2], hi[3])],
                             [(lo[1], hi[2], lo[3]), (lo[1], hi[2], hi[3]), (hi[1], hi[2], hi[3]), (hi[1], hi[2], lo[3])],
                             [(lo[1], lo[2], lo[3]), (lo[1], hi[2], lo[3]), (hi[1], hi[2], lo[3]), (hi[1], lo[2], lo[3])],
                             [(lo[1], lo[2], hi[3]), (hi[1], lo[2], hi[3]), (hi[1], hi[2], hi[3]), (lo[1], hi[2], hi[3])]],
                            [(0, (0, 0, 0)) for _ in 1:6])
    # ball containing the body: exact volume; ball inside: the ball
    @test HR._ball_volume(box((0, 0, 0), (1, 1, 1)), (0.5, 0.5, 0.5), 1.0) ≈ 1.0
    @test HR._ball_volume(box((0, 0, 0), (1, 1, 1)), (0.5, 0.5, 0.5), 0.3) ≈ 4π * 0.3^3 / 3 rtol = 1e-12
    # one plane cutting a cap: sphere minus the cap of height R - h
    R, h = 0.4, 0.15
    v = HR._ball_volume(box((-1, -1, -1), (h, 1, 1)), (0.0, 0.0, 0.0), R)
    a = R - h
    @test v ≈ 4π * R^3 / 3 - π * a^2 * (3R - a) / 3 rtol = 2e-4
    # 2D: disk cut by a line
    sq = HR._Poly2([(-1.0, -1.0), (h, -1.0), (h, 1.0), (-1.0, 1.0)], [(0, (0, 0)) for _ in 1:4])
    seg = R^2 * acos(h / R) - h * sqrt(R^2 - h^2)
    @test HR._ball_volume(sq, (0.0, 0.0), R) ≈ π * R^2 - seg rtol = 1e-5
    # the radius keeps away from tangencies
    @test HR._clear_radius(1.0, [1.001]) ≈ 1.001 * 0.998 rtol = 1e-6
    @test HR._clear_radius(1.0, [0.996, 0.998, 1.0, 1.002, 1.004]) == 1.0      # no room within 0.5 %
    @test HR._clear_radius(1.0, [0.5, 2.0]) == 1.0
end

@testset "rounded pores: density, wall and rounding" begin
    cell = Cell((1.0, 1.0, 1.0))
    t = laguerre(cell; spheres=dense_spheres(cell; radius=lognormal_r, fraction=0.55, rng=3), volumes=:spheres)
    # rounding 0: the polyhedral pores, unchanged
    @test pores(t; wall=0.01, rounding=0) isa Vector{Polyhedron}
    @test wall_thickness(t, 0.12; rounding=0) == wall_thickness(t, 0.12)
    # density increases with rounding; rounding and wall are inverse problems
    ρ(r) = HR._foam_density(t, 0.006, r)
    @test ρ(0.0) <= ρ(0.5) < ρ(0.9) < ρ(1.0)
    r = pore_rounding(t, 0.12; wall=0.006)
    @test 0.9 < r < 1
    @test ρ(r) ≈ 0.12 rtol = 2e-3
    @test wall_thickness(t, 0.12; rounding=r) ≈ 0.006 rtol = 1e-4
    # out of reach
    @test_throws r"below the polyhedral foam" pore_rounding(t, 0.02; wall=0.006)
    @test_throws r"above the fully rounded" pore_rounding(t, 0.5; wall=0.006)
    @test_throws r"reduce rounding" wall_thickness(t, 0.05; rounding=1.0)
    @test_throws r"rounding must be" pores(t; wall=0.01, rounding=1.5)
    @test_throws r"rounding must be" wall_thickness(t, 0.1; rounding=-0.1)
    # rounded pores: balls, intersections or polyhedra; the mesh has the predicted density
    P = pores(t; wall=0.006, rounding=r)
    @test all(p -> p isa Union{Intersect,Polyhedron,Sphere,Rotate}, P)
    @test count(p -> p isa Intersect, P) > length(t.loops) / 2
    res = generate(RVE(cell, Phase(:pores, P; void=true); matrix=:walls); h=0.06, retries=3)
    @test volume_fraction(res, :walls) ≈ 0.12 rtol = 2e-3
    @test all(x -> x.ok, res.periodicity) && res.stats.ninverted == 0
    # 2D
    c2 = Cell((1.0, 1.0))
    t2 = laguerre(c2; spheres=dense_spheres(c2; radius=r -> 0.05 + 0.05rand(r), fraction=0.7, rng=2), volumes=:spheres)
    r2 = pore_rounding(t2, 0.15; wall=0.01)
    res2 = generate(RVE(c2, Phase(:pores, pores(t2; wall=0.01, rounding=r2); void=true); matrix=:walls); h=0.03, retries=3)
    @test volume_fraction(res2, :walls) ≈ 0.15 rtol = 1e-3
    # a block cut from a periodic tessellation (cell not periodic along z)
    cz = Cell((1.0, 1.0, 0.6); periodic=(true, true, false))
    tz = laguerre(cz; seeds=[x.center for x in dense_spheres(Cell((1.0, 1.0, 0.6)); radius=0.12, n=20, rng=4)], periodic=true)
    rz = pore_rounding(tz, 0.2; wall=0.01)
    resz = generate(RVE(cz, Phase(:pores, pores(tz; wall=0.01, rounding=rz); void=true); matrix=:walls); h=0.08, retries=3)
    @test volume_fraction(resz, :walls) ≈ 0.2 rtol = 2e-3
end

@testset "generate retries with a shift" begin
    @test_throws ArgumentError generate(RVE(Cell((1.0, 1.0, 1.0)), Phase(:s, Sphere((0.5, 0.5, 0.5), 0.2))); h=0.3, retries=-1)
    @test_throws ArgumentError generate(RVE(Cell((1.0, 1.0, 1.0)), Phase(:s, Sphere((0.5, 0.5, 0.5), 0.2))); h=0.3, retries=4)
    # a failure on the periodic faces is retried on the shifted material (a stub failure)
    res = generate(RVE(Cell((1.0, 1.0, 1.0)), Phase(:s, Sphere((0.5, 0.5, 0.5), 0.2))); h=0.3, retries=1)
    @test res.stats.nelements > 0
    @test HR._periodic_failure(ErrorException("Gmsh failed to mesh the RVE (Periodic meshing of surface 1 with surface 2)"))
    @test HR._periodic_failure(PeriodicityError(1, "x"))
    @test !HR._periodic_failure(ArgumentError("Periodic"))
end
