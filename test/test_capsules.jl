# Random short fibres (capsules) and orientation tensors.

const SD = HyperRVE._segdist2
# independent check: smallest surface gap between capsules, exact (ternary search of the convex
# distance from a point of one axis to the other segment), all images that can touch
function pseg(p, a, b)
    ab = b .- a; l2 = sum(abs2, ab)
    t = l2 == 0 ? 0.0 : clamp(sum((p .- a) .* ab) / l2, 0, 1)
    return sqrt(sum(abs2, p .- a .- t .* ab))
end
function segseg(a1, b1, a2, b2)
    lo, hi = 0.0, 1.0
    f(t) = pseg(a1 .+ t .* (b1 .- a1), a2, b2)
    for _ in 1:200
        m1 = lo + (hi - lo) / 3; m2 = hi - (hi - lo) / 3
        f(m1) < f(m2) ? (hi = m2) : (lo = m1)
    end
    return min(f(lo), f(0.0), f(1.0))
end
function capsule_min_gap(C, cell; m=nothing)
    L = cell.size; per = cell.periodic; g = Inf
    ext = maximum(maximum(abs.(c.p2 .- c.p1)) + 2c.radius for c in C)
    K = ntuple(d -> per[d] ? ceil(Int, ext / L[d]) + 1 : 0, 3)
    for i in eachindex(C), j in i:length(C), k in Iterators.product(ntuple(d -> -K[d]:K[d], 3)...)
        i == j && all(==(0), k) && continue
        s = k .* L
        g = min(g, segseg(C[i].p1, C[i].p2, C[j].p1 .+ s, C[j].p2 .+ s) - C[i].radius - C[j].radius)
    end
    return g
end

@testset "segment distance" begin
    @test SD((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)) ≈ 1            # parallel
    @test SD((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, -1.0, 2.0), (0.5, 1.0, 2.0)) ≈ 4            # crossing at distance 2
    @test SD((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)) ≈ 1            # collinear
    @test SD((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 1.0, 0.0), (3.0, 5.0, 0.0)) ≈ 2            # end to end
    @test SD((0.5, 3.0, 0.0), (0.5, 3.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)) ≈ 9            # point-segment
    @test SD((1.0, 2.0, 2.0), (1.0, 2.0, 2.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)) ≈ 9            # point-point
    @test SD((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.2, 0.0, 0.0), (0.7, 0.0, 0.0)) ≈ 0            # overlapping
end

@testset "random capsules" begin
    c = unitcell()
    f = random_capsules(c; length=0.3, radius=0.04, fraction=0.03, rng=1)
    @test length(f) == round(Int, 0.03 / (π * 0.04^2 * 0.3 + 4π * 0.04^3 / 3))
    @test all(x -> x.radius == 0.04 && sqrt(sum(abs2, x.p2 .- x.p1)) ≈ 0.3, f)
    @test volume_fraction(f, c) ≈ length(f) * (π * 0.04^2 * 0.3 + 4π * 0.04^3 / 3)
    @test capsule_min_gap(f, c) >= 0.004 - 1e-9
    @test random_capsules(c; length=0.3, radius=0.04, n=8, rng=2) == random_capsules(c; length=0.3, radius=0.04, n=8, rng=2)
    @test random_capsules(c; length=0.3, radius=0.04, n=8, rng=2) != random_capsules(c; length=0.3, radius=0.04, n=8, rng=3)
    # vector of lengths, sampler, sorted
    f = random_capsules(c; length=[0.1, 0.3, 0.2], radius=0.03, rng=1)
    @test [sqrt(sum(abs2, x.p2 .- x.p1)) for x in f] ≈ [0.3, 0.2, 0.1]
    f = random_capsules(c; length=rng -> 0.1 + 0.2rand(rng), radius=0.03, fraction=0.02, rng=4)
    @test 0.015 < volume_fraction(f, c) <= 0.02
    # long fibres in a small, non-cubic cell away from the origin: images of themselves
    cs = Cell((0.5, 0.6, 0.7); origin=(2.0, -1.0, 0.3))
    f = random_capsules(cs; length=0.9, radius=0.03, n=10, rng=5)
    @test length(f) == 10 && capsule_min_gap(f, cs) >= 0.003 - 1e-9
    # mixed periodicity and :inside
    cm = Cell((1.0, 1.0, 0.6); periodic=(true, true, false))
    f = random_capsules(cm; length=0.25, radius=0.03, n=12, rng=6)
    @test capsule_min_gap(f, cm) >= 0.003 - 1e-9
    f = random_capsules(c; length=0.25, radius=0.03, n=12, boundary=:inside, rng=7)
    @test all(x -> all(min.(x.p1, x.p2) .- 0.03 .>= 0.003 - 1e-12) && all(max.(x.p1, x.p2) .+ 0.03 .<= 1 - 0.003 + 1e-12), f)
    # face bands: surface and junction circles away from the face planes
    f = random_capsules(c; length=0.3, radius=0.04, n=15, rng=8)
    for x in f
        u = (x.p2 .- x.p1) ./ 0.3
        for d in 1:3
            e = 0.04 * sqrt(max(0, 1 - u[d]^2))
            for v in (min(x.p1[d], x.p2[d]) - 0.04, max(x.p1[d], x.p2[d]) + 0.04, x.p1[d] - e, x.p1[d] + e, x.p2[d] - e, x.p2[d] + e)
                @test abs(v - round(v)) >= 0.004 - 1e-12
            end
        end
    end
    # avoid spheres (hybrid composite) and avoid by spheres
    p = random_spheres(c; radius=0.12, n=4, rng=1)
    f = random_capsules(c; length=0.3, radius=0.04, n=10, avoid=p, rng=2)
    g = minimum(sqrt(HyperRVE._segdist2(Tuple(x.p1), Tuple(x.p2), s.center .+ k, s.center .+ k)) - x.radius - s.radius
                for x in f for s in p for k in Iterators.product(-1:1, -1:1, -1:1))
    @test g >= 0.004 - 1e-12
    s = random_spheres(c; radius=0.05, n=15, avoid=f, rng=3)
    g = minimum(sqrt(HyperRVE._segdist2(Tuple(x.p1 .+ k), Tuple(x.p2 .+ k), q.center, q.center)) - x.radius - q.radius
                for x in f for q in s for k in Iterators.product(-1:1, -1:1, -1:1))
    @test g >= 0.005 - 1e-12
    # avoid infinite fibres
    fb = random_fibres(c; radius=0.1, n=2, rng=1)
    f = random_capsules(c; length=0.4, radius=0.03, n=10, avoid=fb, rng=4)
    g = minimum(sqrt(sum(abs2, ((x.p1 .+ t .* (x.p2 .- x.p1)) .- y.point .- k)[1:2])) - x.radius - y.radius
                for x in f for y in fb for k in Iterators.product(-1:1, -1:1, 0:0) for t in range(0, 1; length=200))
    @test g >= 0.003 - 0.4 / 199
end

@testset "orientations" begin
    c = Cell((3.0, 3.0, 3.0))
    f = random_capsules(c; length=0.2, radius=0.005, n=2000, rng=1)
    A = orientation_tensor(f)
    @test A ≈ A' && sum(A[i, i] for i in 1:3) ≈ 1
    @test all(abs.(A .- [i == j ? 1 / 3 : 0.0 for i in 1:3, j in 1:3]) .< 0.03)
    f = random_capsules(c; length=0.2, radius=0.005, n=500, orientation=:planar, axis=1, rng=2)
    A = orientation_tensor(f)
    @test A[1, 1] < 1e-12 && abs(A[2, 2] - 0.5) < 0.05 && abs(A[3, 3] - 0.5) < 0.05
    f = random_capsules(c; length=0.2, radius=0.01, n=50, orientation=(1, 1, 0), rng=3)
    @test all(x -> abs(x.p2[3] - x.p1[3]) < 1e-12 && abs((x.p2 .- x.p1)[1] - (x.p2 .- x.p1)[2]) < 1e-12, f)
    @test orientation_tensor(f) ≈ [0.5 0.5 0; 0.5 0.5 0; 0 0 0]
    f = random_capsules(c; length=0.2, radius=0.01, n=20, orientation=rng -> (0, 0, rand(rng) < 0.5 ? 1 : -1), rng=4)
    @test orientation_tensor(f)[3, 3] ≈ 1
    @test orientation_tensor(random_fibres(unitcell(); radius=0.1, n=3, rng=1))[3, 3] ≈ 1
end

@testset "capsule meshes" begin
    c = unitcell()
    f = random_capsules(c; length=0.3, radius=0.05, n=6, rng=11)
    res = generate(RVE(c, Phase(:fibres, f; refine=0.03)); h=0.12)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :fibres) ≈ volume_fraction(f, c) rtol = 2e-3
    @test res.stats.min_quality > 0.01
    f = random_capsules(c; length=0.4, radius=0.05, n=6, orientation=(0, 0, 1), rng=12)
    res = generate(RVE(c, Phase(:fibres, f; refine=0.03)); h=0.12)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :fibres) ≈ volume_fraction(f, c) rtol = 2e-3
    # fibres longer than the cell, and a hybrid with overlaps = :error
    cs = Cell((0.5, 0.6, 0.7))
    f = random_capsules(cs; length=0.9, radius=0.03, n=4, rng=5)
    res = generate(RVE(cs, Phase(:fibres, f; refine=0.02)); h=0.08)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity) && res.stats.min_quality > 0
    p = random_spheres(c; radius=0.12, n=3, rng=4)
    f = random_capsules(c; length=0.3, radius=0.04, n=6, avoid=p, min_gap=0.02, rng=5)
    res = generate(RVE(c, [Phase(:p, p; refine=0.04), Phase(:f, f; refine=0.03)]; overlaps=:error); h=0.12)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity) && res.stats.min_quality > 0
end

@testset "capsule errors" begin
    c = unitcell()
    @test_throws ArgumentError random_capsules(c; length=0.3, n=5)
    @test_throws ArgumentError random_capsules(c; radius=0.03, n=5)
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03)
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, fraction=0.1)
    @test_throws ArgumentError random_capsules(c; length=[0.3], radius=0.03, n=5)
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, orientation=:random)
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, orientation=(0, 0, 0))
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, orientation=[1 0 0; 0 1 0; 0 0 0])
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, orientation=:planar, axis=4)
    @test_throws ArgumentError random_capsules(c; length=1.2, radius=0.03, n=5, orientation=(0, 0, 1))   # use random_fibres
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.6, n=1)
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, boundary=:wall)
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, avoid=Box((0, 0, 0), (0.1, 0.1, 0.1)))
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, avoid=Disk((0.5, 0.5), 0.1))
    @test_throws MethodError random_capsules(Cell((1.0, 1.0)); length=0.3, radius=0.03, n=5)
    @test_throws ErrorException random_capsules(c; length=0.6, radius=0.08, fraction=0.4, rng=1, max_tries=300)
    f = @test_logs (:warn,) random_capsules(c; length=0.6, radius=0.08, fraction=0.4, rng=1, max_tries=300, strict=false)
    @test 0 < length(f)
    @test_throws ArgumentError orientation_tensor([Sphere((0, 0, 0), 1.0)])
    @test_throws ArgumentError orientation_tensor(Capsule[])
end

@testset "capsule review fixes" begin
    # fibres much longer than the cell: never touching their own images, any seed
    c = unitcell()
    for seed in 1:8
        f = random_capsules(c; length=2.3, radius=0.03, n=6, rng=seed, max_tries=3000)
        @test length(f) == 6 && capsule_min_gap(f, c) >= 0.003 - 1e-9
    end
    cs = Cell((0.5, 0.6, 0.7))
    for seed in 1:6
        f = random_capsules(cs; length=0.9, radius=0.03, n=6, rng=seed)
        @test length(f) == 6 && capsule_min_gap(f, cs) >= 0.003 - 1e-9
    end
    # :inside with fibres almost as long as the cell, non-cubic cell away from the origin
    ci = Cell((1.0, 1.2, 0.9); origin=(1.0, 2.0, -3.0))
    for seed in 1:5
        f = random_capsules(ci; length=0.85, radius=0.03, n=5, boundary=:inside, rng=seed)
        @test length(f) == 5
        @test all(x -> all(min.(x.p1, x.p2) .- 0.03 .- ci.origin .>= 0.003 - 1e-12) &&
                       all(max.(x.p1, x.p2) .+ 0.03 .- ci.origin .<= ci.size .- 0.003 .+ 1e-12), f)
    end
    @test_throws ArgumentError random_capsules(c; length=1.8, radius=0.03, n=2, boundary=:inside, rng=1)
    # oblique fixed direction longer than the cell is accepted when it fits
    f = random_capsules(c; length=1.3, radius=0.02, n=3, orientation=(1, 2, 3), rng=1)
    @test length(f) == 3 && capsule_min_gap(f, c) >= 0.002 - 1e-9
    # bands: centres and seams away from the faces, in-plane fibres included
    flat = Cell((1.0, 1.0, 0.4))
    for seed in 1:20
        f = random_capsules(flat; length=0.35, radius=0.04, n=8, orientation=:planar, axis=3, rng=seed)
        for x in f, v in (x.p1[3], x.p2[3], x.p1[3] - 0.04, x.p1[3] + 0.04)
            @test abs(v - 0.4 * round(v / 0.4)) >= 0.004 - 1e-12
        end
    end
    @test_throws ArgumentError random_capsules(c; length=0.3, radius=0.03, n=5, axis=7)
    @test random_capsules(c; length=0.3, radius=0.04, n=5, rng=3) == random_capsules(c; length=0.3, radius=0.04, n=5, rng=HyperRVE.Random.Xoshiro(3))
    # orientation tensor: number average, translates
    @test orientation_tensor([Capsule((0, 0, 0), (2, 0, 0), 0.1), Translate(Capsule((0, 0, 0), (0, 0.1, 0), 0.1), (1, 0, 0))]) ≈
          [0.5 0 0; 0 0.5 0; 0 0 0]
end

