# Particle mixtures (avoid) and clustered packings.

# independent check: smallest surface gap between two sets, all periodic images (3^D)
function cross_gap(A, B, cell)
    D = length(cell.size); L = cell.size; per = cell.periodic
    pt(g) = g isa Fibre ? g.point : g.center
    g = Inf
    for a in A, b in B, k in Iterators.product(ntuple(d -> per[d] ? (-1:1) : (0:0), D)...)
        δ = collect(pt(a) .- pt(b) .+ k .* L)
        for x in (a, b)
            x isa Fibre && (δ[findfirst(v -> abs(v) == 1, x.direction)] = 0)
        end
        g = min(g, sqrt(sum(abs2, δ)) - a.radius - b.radius)
    end
    return g
end
self_gap(S, cell) = minimum(cross_gap([S[i]], [S[j]], cell) for i in eachindex(S) for j in i+1:length(S))
# smallest distance from sphere surfaces to a capsule (dense sampling of the axis), with images
function capsule_gap(S, cap, cell)
    L = cell.size; g = Inf
    for s in S, k in Iterators.product((-1:1, -1:1, -1:1)...), t in range(0, 1; length=2001)
        p = cap.p1 .+ t .* (cap.p2 .- cap.p1) .+ k .* L
        g = min(g, sqrt(sum(abs2, s.center .- p)) - s.radius - cap.radius)
    end
    return g
end

@testset "avoid: populations" begin
    c = unitcell()
    # unchanged stream without avoid / cluster
    @test random_spheres(c; radius=0.08, n=12, avoid=Sphere[], rng=3) == random_spheres(c; radius=0.08, n=12, rng=3)
    @test random_spheres(c; radius=0.08, n=12, avoid=nothing, cluster=nothing, rng=3) == random_spheres(c; radius=0.08, n=12, rng=3)
    hard = random_spheres(c; radius=0.14, fraction=0.10, rng=1)
    soft = random_spheres(c; radius=0.09, fraction=0.06, avoid=hard, rng=2)
    @test volume_fraction(soft, c) ≈ 0.06 rtol = 0.1
    @test cross_gap(hard, soft, c) >= 0.009 - 1e-12                  # min_gap = 0.1 * 0.09
    @test random_spheres(c; radius=0.09, fraction=0.06, avoid=hard, rng=2) == soft
    @test random_spheres(c; radius=0.09, fraction=0.06, avoid=hard, rng=HyperRVE.Random.Xoshiro(2)) == soft
    # a third population avoiding both, as a tuple of collections
    third = random_spheres(c; radius=0.05, n=15, avoid=(hard, soft), min_gap=0.02, rng=5)
    @test cross_gap(hard, third, c) >= 0.02 - 1e-12 && cross_gap(soft, third, c) >= 0.02 - 1e-12
    # cell away from the origin, non-periodic direction
    c2 = Cell((1.0, 0.8, 1.2); origin=(3.0, -1.0, 0.5), periodic=(true, false, true))
    h2 = random_spheres(c2; radius=0.12, n=8, rng=6)
    s2 = random_spheres(c2; radius=0.07, n=20, avoid=h2, rng=7)
    @test cross_gap(h2, s2, c2) >= 0.007 - 1e-12
    # capsules, cylinders (by their capsule), ellipsoids from chains, fibres, translates
    cap = Capsule((0.1, 0.5, 0.5), (0.9, 0.5, 0.5), 0.1)
    s = random_spheres(c; radius=0.06, n=25, avoid=cap, rng=8)
    @test capsule_gap(s, cap, c) >= 0.006 - 2e-4
    cyl = Cylinder((0.5, 0.1, 0.2), (0.0, 0.8, 0.0), 0.1)
    s = random_spheres(c; radius=0.06, n=25, avoid=Translate(cyl, (0.1, 0.0, 0.0)), rng=8)
    @test capsule_gap(s, Capsule((0.6, 0.1, 0.2), (0.6, 0.9, 0.2), 0.1), c) >= 0.006 - 2e-4
    ch = particle_chains(c; radius=0.08, gap=0.02, n_chains=2, aspect=1.5, rng=3)
    @test all(x -> x isa Ellipsoid, ch)
    s = random_spheres(c; radius=0.05, n=20, avoid=ch, rng=9)
    # the ellipsoids lie inside their capsules: check against those
    caps = [let k = argmax(collect(e.radii)), b = maximum(e.radii[i] for i in 1:3 if i != k),
                ev = ntuple(i -> i == k ? e.radii[k] - b : 0.0, 3)
                Capsule(e.center .- ev, e.center .+ ev, b) end for e in ch]
    @test minimum(capsule_gap(s, cp, c) for cp in caps) >= 0.005 - 2e-4
    fb = random_fibres(c; radius=0.1, n=3, rng=1)
    s = random_spheres(c; radius=0.05, n=20, avoid=fb, rng=10)
    @test cross_gap(fb, s, c) >= 0.005 - 1e-12
    # fibres avoiding fibres; 2D disks avoiding disks
    f2 = random_fibres(c; radius=0.05, n=10, avoid=fb, rng=11)
    @test cross_gap(fb, f2, c) >= 0.005 - 1e-12
    d2 = Cell((2.0, 1.0))
    big = random_disks(d2; radius=0.15, n=6, rng=1)
    small = random_disks(d2; radius=0.05, fraction=0.1, avoid=big, rng=2)
    @test cross_gap(big, small, d2) >= 0.005 - 1e-12
end

@testset "clusters" begin
    c = unitcell()
    c0 = (0.5, 0.5, 0.5)
    s = random_spheres(c; radius=0.04, n=12, cluster=(centers=[c0], radius=0.25), rng=1)
    mind(x) = sqrt(sum(abs2, (x .- c0) .- round.(x .- c0)))
    @test length(s) == 12 && all(x -> mind(x.center) <= 0.25 + 1e-12, s)
    # clustered packings have closer neighbours than uniform ones at the same fraction
    nn(S) = sum(minimum(sqrt(sum(abs2, (a.center .- b.center) .- round.(a.center .- b.center))) for b in S if b !== a) for a in S) / length(S)
    u = random_spheres(c; radius=0.04, n=40, rng=2)
    cl = random_spheres(c; radius=0.04, n=40, cluster=(n=3, radius=0.2), rng=2)
    @test nn(cl) < nn(u)
    @test self_gap(cl, c) >= 0.004 - 1e-12
    @test random_spheres(c; radius=0.04, n=40, cluster=(n=3, radius=0.2), rng=2) == cl
    # share < 1 puts some inclusions outside the balls
    sh = random_spheres(c; radius=0.03, n=40, cluster=(centers=[c0], radius=0.15, share=0.5), rng=3)
    @test any(x -> mind(x.center) > 0.15 + 1e-9, sh) && any(x -> mind(x.center) <= 0.15, sh)
    # 2D, cluster across the periodic faces
    d2 = Cell((1.0, 1.0); origin=(-2.0, 1.0))
    dk = random_disks(d2; radius=0.03, n=15, cluster=(centers=[(-2.0, 1.0)], radius=0.2), rng=4)
    @test all(x -> sqrt(sum(abs2, (x.center .- (-2.0, 1.0)) .- round.(x.center .- (-2.0, 1.0)))) <= 0.2 + 1e-12, dk)
    @test self_gap(dk, d2) >= 0.003 - 1e-12
end

@testset "mixture meshes" begin
    c = unitcell()
    hard = random_spheres(c; radius=0.14, n=3, rng=1)
    soft = random_spheres(c; radius=0.08, n=6, avoid=hard, rng=2)
    res = generate(RVE(c, [Phase(:hard, hard), Phase(:soft, soft)]; overlaps=:error); h=0.12)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :hard) ≈ volume_fraction(hard, c) rtol = 1e-3
    @test volume_fraction(res, :soft) ≈ volume_fraction(soft, c) rtol = 1e-3
    @test res.stats.min_quality > 0
    @test !("interface_hard_soft" in interface_names(res))
    d2 = Cell((1.0, 1.0))
    cl = random_disks(d2; radius=0.04, fraction=0.12, cluster=(n=3, radius=0.2), rng=5)
    un = random_disks(d2; radius=0.025, fraction=0.05, avoid=cl, rng=6)
    res = generate(RVE(d2, [Phase(:cl, cl), Phase(:un, un)]; overlaps=:error); h=0.04)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :un) ≈ volume_fraction(un, d2) rtol = 1e-3
end

@testset "avoid and cluster errors" begin
    c = unitcell()
    @test_throws ArgumentError random_spheres(c; radius=0.05, n=5, avoid=Box((0, 0, 0), (0.2, 0.2, 0.2)))
    @test_throws ArgumentError random_spheres(c; radius=0.05, n=5, avoid=Fibre((0.5, 0.5, 0.5), (1, 1, 0), 0.1))
    @test_throws ArgumentError random_spheres(c; radius=0.05, n=5, avoid=Disk((0.5, 0.5), 0.1))
    @test_throws ArgumentError random_spheres(c; radius=0.05, n=5, avoid=[1, 2])
    @test_throws ArgumentError random_fibres(c; radius=0.05, n=5, avoid=Sphere((0.5, 0.5, 0.5), 0.1))
    @test_throws ArgumentError random_fibres(c; radius=0.05, n=5, avoid=Fibre((0.5, 0.5, 0.5), (1, 0, 0), 0.1))
    for cl in ((n=2,), (radius=0.2,), (n=2, centers=[(0.5, 0.5, 0.5)], radius=0.2), (n=0, radius=0.2), (n=2, radius=-1),
               (n=2, radius=0.2, share=0), (n=2, radius=0.2, share=1.5), (n=2, radius=0.2, spread=1), (2, 0.2))
        @test_throws ArgumentError random_spheres(c; radius=0.05, n=5, cluster=cl, rng=1)
    end
    # too many inclusions for the clusters
    @test_throws ArgumentError random_spheres(c; radius=0.05, n=60, cluster=(n=1, radius=0.1), rng=1)
    # saturation around obstacles: error, or warning and partial packing
    big = Sphere((0.5, 0.5, 0.5), 0.45)
    @test_throws ErrorException random_spheres(c; radius=0.1, n=40, avoid=big, rng=1, max_tries=2000)
    s = @test_logs (:warn,) random_spheres(c; radius=0.1, n=40, avoid=big, rng=1, max_tries=2000, strict=false)
    @test 0 < length(s) < 40
end
