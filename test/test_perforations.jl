using HyperRVE, Test, Random

@testset "perforate" begin
    dist = truncated_lognormal(4.0; cv=0.35, min=2.0, max=8.0)
    c = Cell((10.0, 10.0, 10.0))
    f = laguerre_foam(c; diameter=dist, porosity=0.9, rng=1)                 # polyhedral pores
    fr = laguerre_foam(c; diameter=dist, porosity=0.9, wall=0.1, rng=1)      # rounded pores
    N = length(pores(f))
    # one hole per cell, centred, along z
    a = perforate(f; diameter=1.5)
    r = a.report
    @test r.holes == 1 == length(a.holes) && r.spacing == (10.0, 10.0) && r.axis == 3
    @test only(a.holes).direction == (0.0, 0.0, 1.0) && only(a.holes).radius == 0.75
    @test only(r.positions) == (5.0, 5.0, 0.0)
    @test r.perforation_ratio ≈ π * 0.75^2 / 100
    @test r.open_pores == r.crossed_pores > 0 && r.open_pores + r.closed_pores == N
    @test r.networks == 1 && r.percolates == (false, false, true) && r.windows == 0
    @test r.porosity ≈ 0.9 rtol = 1e-6
    @test 0 < r.open_porosity < r.porosity
    @test length(a.fluid) + length(a.closed) == N + 1
    # exact volumes (flat walls): solid + air + closed pores = cell, closed = porosity - open porosity
    g = generate(RVE(c, foam_phases(a); matrix=:walls); geometry_only=true)
    @test 1 - volume_fraction(g, :air) - volume_fraction(g, :walls) ≈ r.porosity - r.open_porosity rtol = 1e-6
    @test volume_fraction(g, :air) > r.open_porosity
    # 2×2 holes; more pores open; deterministic
    a4 = perforate(f; diameter=1.0, spacing=5.0)
    @test a4.report.holes == 4 && a4.report.positions[1] == (2.5, 2.5, 0.0)
    @test a4.report.open_pores >= r.open_pores
    @test perforate(f; diameter=1.0, spacing=(5.0, 10.0)).report.holes == 2
    # windows join the networks of separate holes
    if a4.report.networks > 1
        @test perforate(f; diameter=1.0, spacing=5.0, windows=(open_fraction=0.4,), rng=3).report.networks < a4.report.networks
    end
    aw = perforate(f; diameter=1.5, windows=(open_fraction=0.3,), rng=3)
    @test aw.report.windows > 0 && aw.report.open_pores >= r.open_pores
    @test aw.report.percolates[3] && aw.report.networks == 1
    @test length(aw.windows) == aw.report.windows == length(aw.report.window_diameters)
    # isolated pores kept on request
    ak = perforate(f; diameter=1.5, isolated=:keep)
    @test isempty(ak.closed)
    @test ak.report.open_porosity ≈ 0.9 rtol = 1e-6
    # other axes, random and given positions
    ax = perforate(f; diameter=1.5, axis=1)
    @test only(ax.holes).direction == (1.0, 0.0, 0.0) && ax.report.percolates == (true, false, false)
    @test only(ax.report.positions) == (0.0, 5.0, 5.0)
    pr = perforate(f; diameter=1.2, position=:random, rng=2)
    @test pr.report.positions == perforate(f; diameter=1.2, position=:random, rng=2).report.positions
    for p in (only(pr.report.positions)[1], only(pr.report.positions)[2]), x in (p, 10 - p)
        @test !(0.3 < x < 0.75)
    end
    p0 = perforate(f; diameter=1.5, position=(0.0, 0.0))             # on the cell edge
    @test only(p0.report.positions) == (0.0, 0.0, 0.0) && p0.report.crossed_pores > 0
    g0 = generate(RVE(c, foam_phases(p0); matrix=:walls); geometry_only=true)
    @test 1 - volume_fraction(g0, :air) - volume_fraction(g0, :walls) ≈ p0.report.porosity - p0.report.open_porosity rtol = 1e-6
    # origin and units do not matter
    k = 1e-3
    ck = Cell(k .* c.size; origin=(1.0, 2.0, 3.0))
    fk = Foam{3}(ck, laguerre(ck; seeds=[ck.origin .+ k .* x .* c.size for x in f.tessellation.seeds],
                 volumes=k^3 .* grain_volumes(f.tessellation)), f.porosity, k * f.wall, f.rounding, k .* f.diameters, f.scale)
    rk = perforate(fk; diameter=k * 1.0, spacing=k * 5.0).report
    @test (rk.crossed_pores, rk.open_pores, rk.networks) == (a4.report.crossed_pores, a4.report.open_pores, a4.report.networks)
    @test collect(rk.positions[1]) ≈ collect((1.0, 2.0, 3.0) .+ k .* (2.5, 2.5, 0.0))
    # a cell not periodic along the holes: through-holes between the free faces
    cz = Cell((10.0, 10.0, 10.0); periodic=(true, true, false))
    fz = laguerre_foam(cz; diameter=dist, porosity=0.9, rng=1)
    az = perforate(fz; diameter=1.5)
    @test az.report.percolates == (false, false, false) && az.report.open_pores > 0
    gz = generate(RVE(cz, foam_phases(az); matrix=:walls); geometry_only=true)
    @test 1 - volume_fraction(gz, :air) - volume_fraction(gz, :walls) ≈ az.report.porosity - az.report.open_porosity rtol = 1e-6
    # rounded pores, with windows: the mesh is periodic
    ar = perforate(fr; diameter=1.5, windows=(open_fraction=0.2,), rng=3)
    res = generate(RVE(c, foam_phases(ar); matrix=:aluminium); h=1.5, retries=3)
    @test all(x -> x.ok, res.periodicity)
    @test volume_fraction(res, :air) > ar.report.open_porosity
    @test volume_fraction(res, :air) + volume_fraction(res, :aluminium) ≈ 1 - (ar.report.porosity - ar.report.open_porosity) rtol = 5e-3
    # foam_phases also for open_pores
    @test length(foam_phases(ak)) == 1 && length(foam_phases(a)) == 2
    @test foam_phases(open_pores(f; open_fraction=0.3, rng=3); air=:gas)[1].name == "gas"
    # errors
    @test_throws r"does not divide.*3.33333 or 2.5" perforate(f; diameter=1.0, spacing=3.0)
    @test_throws r"too large" perforate(f; diameter=5.0, spacing=5.0)
    @test_throws r"does not mesh" perforate(f; diameter=1.5, position=(0.5, 5.0))
    @test_throws r"position" perforate(f; diameter=1.5, position=:edge)
    @test_throws r"axis" perforate(f; diameter=1.5, axis=4)
    @test_throws r"diameter" perforate(f; diameter=-1.0)
    @test_throws r"isolated" perforate(f; diameter=1.5, isolated=:drop)
    @test_throws r"unknown windows" perforate(f; diameter=1.5, windows=(fraction=0.2,))
    @test_throws r"named tuple" perforate(f; diameter=1.5, windows=0.2)
    @test_throws r"open_fraction" perforate(f; diameter=1.5, windows=(open_fraction=2.0,))
    @test_throws r"only 3D" perforate(laguerre_foam(Cell((20.0, 20.0)); diameter=3.0, porosity=0.85, rng=1); diameter=1.0)
    @test_logs (:warn, r"wider than the mean pore") match_mode = :any perforate(f; diameter=4.5, spacing=10.0)
    @test_throws r"fluid is empty" foam_phases((fluid=Geometry{3}[], closed=Geometry{3}[]))
    @test_throws r"diameter" perforate(f; diameter=true)
    # spacing snapped to divide the cell exactly; holes on a free face: half of their section
    a9 = perforate(f; diameter=1.0, spacing=3.3333333)
    @test a9.report.holes == 9 && a9.report.spacing == (10 / 3, 10 / 3)
    @test a9.report.positions[end][1] ≈ 2.5 * (10 / 3)
    cx = Cell((10.0, 10.0, 10.0); periodic=(false, true, true))
    fx = laguerre_foam(cx; diameter=dist, porosity=0.9, rng=1)
    @test perforate(fx; diameter=1.5, position=(0.0, 5.0)).report.perforation_ratio ≈ π * 0.75^2 / 2 / 100 rtol = 1e-4
    @test only(perforate(f; diameter=1.5, position=(-1e-17, 5.0)).report.positions)[1] == 0.0
    # axis 2 in a non-cubic cell: spacing in the order of the in-plane directions (1, 3)
    cb = Cell((12.0, 8.0, 6.0))
    fb = laguerre_foam(cb; diameter=dist, porosity=0.9, rng=2)
    ab2 = perforate(fb; diameter=1.0, spacing=(4.0, 3.0), axis=2)
    @test ab2.report.holes == 6 && ab2.report.positions[1] == (2.0, 0.0, 1.5)
    @test all(h -> h.direction == (0.0, 1.0, 0.0), ab2.holes) && ab2.report.percolates[2]
    # :keep with windows: every pore in the fluid, networks counts those connected to holes
    awk = perforate(f; diameter=1.5, windows=(open_fraction=0.05,), isolated=:keep, rng=3)
    @test isempty(awk.closed) && awk.report.open_pores == N && awk.report.networks >= 1
    # hole-window edges: a hole image across x = L touches a window of grain 1 (offsets)
    adj = [Tuple{Int,NTuple{3,Int}}[] for _ in 1:N+1]
    win = (1, 2, (0, 1, 0), (9.9, 5.0, 5.0), (0.0, 0.0, 1.0), 0.3, 0.2, 0.2)
    HyperRVE._hole_window_edges!(adj, f.tessellation, [win], [(0.2, 5.0, 0.0)], 0.5, 3)
    @test Set(adj[N+1]) == Set([(1, (-1, 0, 0)), (2, (-1, 1, 0))])
    @test (N + 1, (1, 0, 0)) in adj[1] && (N + 1, (1, -1, 0)) in adj[2]
    adj0 = [Tuple{Int,NTuple{3,Int}}[] for _ in 1:N+1]
    HyperRVE._hole_window_edges!(adj0, f.tessellation, [win], [(5.0, 5.0, 0.0)], 0.5, 3)
    @test isempty(adj0[N+1])
    # holes much wider than the pores: small pores inside a hole are crossed (not left closed over
    # the air): the air with the solid equals the air alone
    fs = laguerre_foam(c; diameter=truncated_lognormal(1.6; cv=0.6, min=0.5, max=6.0), porosity=0.85, rng=3)
    as = @test_logs (:warn, r"wider") match_mode = :any perforate(fs; diameter=4.0, position=(6.984, 4.963))
    v1 = volume_fraction(generate(RVE(c, foam_phases(as); matrix=:walls); geometry_only=true), :air)
    v2 = volume_fraction(generate(RVE(c, Phase(:air, as.fluid); matrix=nothing); geometry_only=true), :air)
    @test abs(v1 - v2) < 1e-4          # before: 1.8e-4 (two pores inside the hole); grazed lenses ~4e-5
    # retries: a hole on the cell edge shifted into the band that does not mesh is skipped
    HyperRVE._FORCED_FAILURES[] = 1
    pe = perforate(f; diameter=0.64, position=(0.0, 0.0))
    re = @test_logs (:warn, r"shifted") match_mode = :any generate(RVE(c, foam_phases(pe); matrix=:walls); h=1.0, retries=3)
    HyperRVE._FORCED_FAILURES[] = 0
    @test all(x -> x.ok, re.periodicity)
end
