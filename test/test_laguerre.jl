using HyperRVE, Test, Random
const HL = HyperRVE

# smallest surface distance between spheres, periodic images along the periodic directions
function min_clearance(s, cell)
    L = cell.size; per = cell.periodic; D = length(L)
    m = Inf
    for i in eachindex(s), j in i+1:length(s)
        δ = ntuple(d -> (x = s[i].center[d] - s[j].center[d]; per[d] ? x - L[d] * round(x / L[d]) : x), D)
        m = min(m, sqrt(sum(abs2, δ)) - s[i].radius - s[j].radius)
    end
    return m
end
# the same with all periodic images (spheres large compared with the cell)
function min_clearance_images(s, cell)
    L = cell.size; per = cell.periodic; D = length(L)
    m = Inf
    for i in eachindex(s), j in i:length(s), k in Iterators.product(ntuple(d -> per[d] ? (-2:2) : (0:0), D)...)
        i == j && all(==(0), k) && continue
        δ = ntuple(d -> s[i].center[d] - s[j].center[d] - k[d] * L[d], D)
        m = min(m, sqrt(sum(abs2, δ)) - s[i].radius - s[j].radius)
    end
    return m
end
# truncated log-normal radius sampler
lognormal(m, cv, lo, hi) = rng -> (s2 = log(1 + cv^2); while true
    d = exp(log(m) - s2 / 2 + sqrt(s2) * randn(rng)); lo <= d <= hi && return d / 2
end)

@testset "dense_spheres" begin
    c = Cell((10.0, 10.0, 10.0))
    s = dense_spheres(c; radius=lognormal(2.0, 0.35, 1.0, 4.0), fraction=0.55, rng=1)
    @test s isa Vector{Sphere}
    @test 0.53 < volume_fraction(s, c) <= 0.55
    @test min_clearance(s, c) >= 0
    @test all(x -> 0.5 <= x.radius <= 2.0, s)                # the radii are not rescaled
    # same rng, same packing
    @test [x.center for x in dense_spheres(c; radius=0.9, n=60, rng=3)] == [x.center for x in dense_spheres(c; radius=0.9, n=60, rng=3)]
    # min_gap, 2D, non-periodic direction (spheres stay inside), origin
    s = dense_spheres(c; radius=0.8, fraction=0.45, min_gap=0.1, rng=2)
    @test min_clearance(s, c) >= 0.1 - 1e-9
    c2 = Cell((1.0, 1.0); origin=(2.0, -1.0))
    d = dense_spheres(c2; radius=0.03, fraction=0.75, rng=1)
    @test d isa Vector{Disk} && volume_fraction(d, c2) > 0.72 && min_clearance(d, c2) >= 0
    cz = Cell((5.0, 5.0, 8.0); periodic=(true, true, false))
    s = dense_spheres(cz; radius=0.5, fraction=0.5, rng=2)
    @test all(x -> x.radius - 1e-9 <= x.center[3] - cz.origin[3] <= 8 - x.radius + 1e-9, s)
    @test min_clearance(s, cz) >= 0
    # beyond jamming: an error, or a warning and smaller spheres
    @test_throws r"jammed at fraction" dense_spheres(Cell((4.0, 4.0, 4.0)); radius=0.5, fraction=0.72, rng=1)
    s = @test_logs (:warn, r"jammed") dense_spheres(Cell((4.0, 4.0, 4.0)); radius=0.5, fraction=0.72, rng=1, strict=false)
    @test 0.55 < volume_fraction(s, Cell((4.0, 4.0, 4.0))) < 0.72 && min_clearance(s, Cell((4.0, 4.0, 4.0))) >= 0
    # min_gap kept when jammed; a gap the smallest spheres cannot keep is an error
    cg = Cell((5.0, 5.0, 5.0))
    s = @test_logs (:warn, r"jammed") dense_spheres(cg; radius=[fill(1.0, 15); fill(0.3, 30)], min_gap=0.15, rng=1, strict=false)
    @test min_clearance(s, cg) >= 0.15 - 1e-9 && all(x -> x.radius > 0, s)
    @test_throws r"min_gap is larger" dense_spheres(cg; radius=[fill(1.0, 15); fill(0.02, 30)], min_gap=0.4, rng=1, strict=false)
    # spheres large compared with the cell: no overlap with any image
    s = dense_spheres(Cell((1.0, 1.0, 1.0)); radius=0.3, n=3, rng=1)
    @test min_clearance_images(s, Cell((1.0, 1.0, 1.0))) >= 0
    @test_throws ArgumentError dense_spheres(c; radius=6.0, n=1)
    @test_throws ArgumentError dense_spheres(c; radius=0.5, fraction=0.3, min_gap=-1)
    # the relaxation does not allocate in its sweeps
    X = [ntuple(_ -> 10rand(), 3) for _ in 1:30]; R = fill(1.0, 30)
    HL._relax!(X, R, (10.0, 10.0, 10.0), (true, true, true); tol=1e-12, maxsweeps=2)
    @test @allocated(HL._relax!(X, R, (10.0, 10.0, 10.0), (true, true, true); tol=1e-12, maxsweeps=20)) < 4096
end

@testset "laguerre tessellations" begin
    c = Cell((10.0, 10.0, 10.0))
    s = dense_spheres(c; radius=lognormal(2.0, 0.35, 1.0, 4.0), fraction=0.55, rng=1)
    R = [x.radius for x in s]
    t = laguerre(c; spheres=s)
    @test t isa Tessellation{3} && t.weights ≈ R .^ 2
    @test sum(grain_volumes(t)) ≈ 1000 rtol = 1e-12
    # volumes = :spheres: exact, proportional to r³
    tv = laguerre(c; spheres=s, volumes=:spheres)
    @test grain_volumes(tv) ≈ R .^ 3 ./ sum(R .^ 3) .* 1000 rtol = 1e-6
    @test grain_diameters(tv) ≈ (6 .* grain_volumes(tv) ./ π) .^ (1 / 3)
    # big spheres get big grains even without the fit
    v = grain_volumes(t)
    @test sum(v[R .> sort(R)[end ÷ 2]]) > sum(v[R .<= sort(R)[end ÷ 2]])
    # equal radii = Voronoi of the same seeds
    pts = [x.center for x in s[1:20]]
    @test grain_volumes(laguerre(c; seeds=pts, radii=0.7)) ≈ grain_volumes(voronoi(c; seeds=pts)) rtol = 1e-9
    @test grain_volumes(laguerre(c; seeds=pts)) ≈ grain_volumes(voronoi(c; seeds=pts)) rtol = 1e-9
    # given volumes
    vt = collect(range(1.0, 3.0; length=20)); vt .*= 1000 / sum(vt)
    @test grain_volumes(laguerre(c; seeds=pts, volumes=vt)) ≈ vt rtol = 1e-6
    @test_throws r"add up" laguerre(c; seeds=pts, volumes=fill(1.0, 20))
    @test_throws r"positive radii" laguerre(c; seeds=pts, volumes=:spheres)
    @test_throws ArgumentError laguerre(c; spheres=s, seeds=pts)
    @test_throws r"one radius per seed" laguerre(c; seeds=pts, radii=[1.0, 2.0])
    # volumes that add up to the cell only to 1e-8; initial weights with an empty cell
    @test grain_volumes(laguerre(c; seeds=pts, volumes=vt .* (1 + 1e-8))) ≈ vt rtol = 1e-6
    c4 = Cell((4.0, 4.0, 4.0))
    sd = [(2.0, 2.0, 2.0), (2.3, 2.0, 2.0), (0.5, 0.5, 0.5)]
    @test grain_volumes(laguerre(c4; seeds=sd, radii=[1.9, 0.0, 0.5], volumes=[20.0, 20.0, 24.0])) ≈ [20, 20, 24] rtol = 1e-6
    # any scale and origin
    for k in (1e-6, 1e4)
        ck = Cell(k .* (10.0, 10.0, 10.0); origin=(3k, -k, 0.0))
        sk = [Sphere(ck.origin .+ k .* (x.center .- c.origin), k * x.radius) for x in s]
        tk = laguerre(ck; spheres=sk, volumes=:spheres)
        @test grain_volumes(tk) ≈ k^3 .* grain_volumes(tv) rtol = 1e-6
        @test wall_thickness(tk, 0.1) ≈ k * wall_thickness(tv, 0.1) rtol = 1e-6
    end
    # a dominated seed: empty grain
    @test_throws r"is empty" laguerre(c4; seeds=[(2.0, 2.0, 2.0), (2.3, 2.0, 2.0)], radii=[1.9, 0.0])
    # non-periodic direction
    cz = Cell((5.0, 5.0, 8.0); periodic=(true, true, false))
    sz = dense_spheres(cz; radius=0.5, fraction=0.5, rng=2)
    tz = laguerre(cz; spheres=sz, volumes=:spheres)
    @test sum(grain_volumes(tz)) ≈ 200 rtol = 1e-12
    @test grain_volumes(tz) ≈ [x.radius^3 for x in sz] ./ sum(x.radius^3 for x in sz) .* 200 rtol = 1e-6
    # 2D
    c2 = Cell((1.0, 1.0))
    d = dense_spheres(c2; radius=r -> 0.02 + 0.03rand(r), n=40, rng=4)
    t2 = laguerre(c2; spheres=d, volumes=:spheres)
    @test grain_volumes(t2) ≈ [x.radius^2 for x in d] ./ sum(x.radius^2 for x in d) rtol = 1e-6
    @test grain_diameters(t2) ≈ sqrt.(4 .* grain_volumes(t2) ./ π)
end

@testset "laguerre foams and meshes" begin
    # closed cells: the walls are oriented by the grain centroid (a Laguerre seed may lie
    # outside its grain), so the porosity of wall_thickness is exact
    c = Cell((10.0, 10.0, 10.0))
    s = dense_spheres(c; radius=lognormal(2.0, 0.35, 1.0, 4.0), fraction=0.55, rng=4)
    t = laguerre(c; spheres=s, volumes=:spheres)
    outside = count(eachindex(s)) do i
        P, planes = HL._grain_planes(t, i)
        x = (s[i].center .- c.origin)
        any(((n, cc),) -> sum(n .* x) > cc + 1e-9, planes)
    end
    @info "Laguerre seeds outside their grain: $outside of $(length(s))"
    w = wall_thickness(t, 0.1)
    @test sum(HL._poly_volume(HL._shrunk(t, i, w / 2)) for i in eachindex(s)) ≈ 900 rtol = 1e-8
    # a seed outside its own grain (a small radius next to a big one)
    c4 = Cell((4.0, 4.0, 4.0))
    t4 = laguerre(c4; seeds=[(1.0, 2.0, 2.0), (1.6, 2.0, 2.0), (3.0, 0.5, 0.7)], radii=[1.0, 0.0, 0.6])
    P, planes = HL._grain_planes(t4, 2)
    @test sum(planes[1][1] .* (1.6, 2.0, 2.0)) > planes[1][2] || any(((n, cc),) -> sum(n .* (1.6, 2.0, 2.0)) > cc + 1e-9, planes)
    w4 = wall_thickness(t4, 0.2)
    @test sum(HL._poly_volume(HL._shrunk(t4, i, w4 / 2)) for i in 1:3) ≈ 0.8 * 64 rtol = 1e-8
    # a periodic Laguerre polycrystal meshes
    cm = Cell((1.0, 1.0, 1.0))
    sm = dense_spheres(cm; radius=r -> 0.12 + 0.1rand(r), n=10, rng=5)
    tm = laguerre(cm; spheres=sm, volumes=:spheres)
    res = generate(RVE(cm, grain_phases(tm); matrix=nothing); h=0.15)
    @test all(x -> x.ok, res.periodicity)
    @test [p.fraction for p in res.phases] ≈ grain_volumes(tm) rtol = 1e-6
end
