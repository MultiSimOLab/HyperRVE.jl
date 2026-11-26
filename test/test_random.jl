# Random packings: keyword combinations, non-overlap with gaps, face bands, reproducibility,
# errors, meshing.

function min_image_gap(items, cell)
    L = cell.size; per = cell.periodic; g = Inf
    c(i) = i isa Fibre ? i.point : i.center
    for a in 1:length(items), b in a+1:length(items)
        δ = collect(c(items[a]) .- c(items[b]))
        for d in eachindex(δ)
            per[d] && (δ[d] -= L[d] * round(δ[d] / L[d]))
            items[a] isa Fibre && items[a].direction[d] == 1 && (δ[d] = 0)
        end
        g = min(g, sqrt(sum(abs2, δ)) - items[a].radius - items[b].radius)
    end
    return g
end

@testset "keyword combinations (R1)" begin
    c = unitcell()
    s = random_spheres(c; radius=0.08, fraction=0.15, rng=1)
    @test length(s) == round(Int, 0.15 / (4π * 0.08^3 / 3)) && all(x -> x.radius == 0.08, s)
    s = random_spheres(c; fraction=0.15, n=20, rng=1)
    @test length(s) == 20 && volume_fraction(s, c) ≈ 0.15
    s = random_spheres(c; radius=[0.1, 0.12, 0.05], rng=1)
    @test [x.radius for x in s] == [0.12, 0.1, 0.05]
    s = random_spheres(c; radius=rng -> 0.05 + 0.03rand(rng), fraction=0.1, rng=2)
    @test 0.08 < volume_fraction(s, c) <= 0.1 && issorted([x.radius for x in s]; rev=true)
    @test random_spheres(c; radius=0.08, n=10, rng=3) == random_spheres(c; radius=0.08, n=10, rng=3)
    @test random_spheres(c; radius=0.08, n=10, rng=3) != random_spheres(c; radius=0.08, n=10, rng=4)
    @test_throws ArgumentError random_spheres(c; radius=0.1)
    @test_throws ArgumentError random_spheres(c; radius=0.1, fraction=0.1, n=3)
    @test_throws ArgumentError random_spheres(c; radius=[0.1], n=3)
    @test_throws ArgumentError random_spheres(c; radius=0.1, fraction=1.2)
    @test_throws ArgumentError random_spheres(c; radius=0.6, n=1)
    @test_throws ArgumentError random_spheres(c; radius=0.1, n=2, boundary=:open)
    @test_throws ArgumentError random_spheres(c; radius=0.05, fraction=0.45)          # above saturation, before trying
    @test_throws ArgumentError random_spheres(c; radius=0.1, n=2, rng="seed")
end

@testset "gaps, bands and boundaries (R2, R3)" begin
    c = Cell((1.0, 2.0, 1.5))
    s = random_spheres(c; radius=0.1, fraction=0.2, rng=5)
    @test min_image_gap(s, c) >= 0.01 - 1e-12
    # face, edge and corner bands, measured on the distances to the planes
    for x in s
        q = x.center .- c.origin
        t = min.(q, c.size .- q)
        @test all(abs.(t .- x.radius) .>= 0.01 - 1e-12) && all(t .>= 0.01 - 1e-12)
        @test abs(sqrt(sum(abs2, t)) - x.radius) >= 0.01 - 1e-12
    end
    s = random_spheres(c; radius=0.1, fraction=0.15, boundary=:inside, rng=6)
    @test all(x -> all(x.center .- x.radius .>= c.origin .+ (0.01 - 1e-12)) &&
                   all(x.center .+ x.radius .<= c.origin .+ c.size .- (0.01 - 1e-12)), s)
    # mixed periodicity: no wrap-around along the non-periodic direction
    cm = Cell((1.0, 1.0, 1.0); periodic=(true, true, false))
    s = random_spheres(cm; radius=0.1, fraction=0.2, min_gap=0.02, rng=7)
    @test min_image_gap(s, cm) >= 0.02 - 1e-12
    d = random_disks(Cell((2.0, 1.0)); radius=0.05, fraction=0.4, rng=8)
    @test min_image_gap(d, Cell((2.0, 1.0))) >= 0.005 - 1e-12 && volume_fraction(d, Cell((2.0, 1.0))) > 0.39
    f = random_fibres(Cell((1.0, 1.0, 3.0)); axis=2, radius=0.06, fraction=0.3, rng=9)
    @test all(x -> x.direction == (0.0, 1.0, 0.0) && x.point[2] == 0.5, f)
    @test volume_fraction(f, Cell((1.0, 1.0, 3.0))) ≈ length(f) * π * 0.06^2 / 3
    @test min_image_gap(f, Cell((1.0, 1.0, 3.0))) >= 0.006 - 1e-12
end

@testset "saturation errors (R4)" begin
    c = unitcell()
    @test_throws ErrorException random_spheres(c; radius=0.08, fraction=0.33, rng=1, max_tries=2000)
    s = @test_logs (:warn, r"placed") match_mode=:any random_spheres(c; radius=0.08, fraction=0.33, rng=1, max_tries=2000, strict=false)
    @test 0.2 < volume_fraction(s, c) < 0.33
end

@testset "meshed packings (R5)" begin
    c = unitcell()
    res = generate(RVE(c, Phase(:p, random_spheres(c; radius=0.12, n=12, rng=11))); h=0.08)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :p) ≈ 12 * 4π * 0.12^3 / 3 rtol = 1e-3
    c2 = Cell((1.0, 1.0))
    res = generate(RVE(c2, Phase(:d, random_disks(c2; radius=0.07, fraction=0.35, rng=12))); h=0.04)
    @test all(x -> x.ok, res.periodicity)
end

struct _ConstRadius end
(::_ConstRadius)(rng) = 0.05

@testset "review fixes" begin
    c = unitcell()
    # volume_fraction: wrong dimensions, oblique fibres, mixed vectors
    @test_throws ArgumentError volume_fraction([Fibre((0.5, 0.5, 0.5), (1, 1, 0), 0.1)], c)
    @test_throws ArgumentError volume_fraction([Disk((0.5, 0.5), 0.1)], c)
    @test_throws ArgumentError volume_fraction([Sphere((0.5, 0.5, 0.5), 0.1)], Cell((1.0, 1.0)))
    @test volume_fraction(Geometry{3}[Sphere((0.2, 0.2, 0.2), 0.1), Fibre((0.6, 0.6, 0.5), (0, 0, 1), 0.1)], c) ≈ 4π * 0.001 / 3 + π * 0.01
    # few inclusions are not rejected by the asymptotic saturation check
    @test length(random_spheres(Cell((1.0, 1.0, 1.0); periodic=false); radius=0.1, fraction=0.4, min_gap=0, face_gap=0, rng=1)) == 95
    @test length(random_spheres(c; fraction=0.4, n=1, rng=1)) == 1
    @test length(random_spheres(c; radius=0.37, n=2, rng=1, face_gap=0)) == 2
    # malformed keywords are ArgumentErrors
    for kw in ((radius=0.1, n=2.5), (radius=0.1, fraction="a"), (radius=0.1, n=2, min_gap=:foo),
               (radius=0.1, n=2, max_tries=0), (radius=0.1, n=2, max_tries=1.5), (radius=0.1, n=2, face_gap=Inf),
               (radius=0.3, n=1, boundary=:inside, face_gap=0.3))
        @test_throws ArgumentError random_spheres(c; rng=1, kw...)
        nothing
    end
    @test_throws ArgumentError random_spheres(c; radius=0.1, n=2, rng=-1)
    @test_throws ArgumentError random_fibres(c; axis=true, radius=0.1, n=2)
    @test_throws ArgumentError random_fibres(c; axis=3.0, radius=0.1, n=2)
    # callable samplers that are not Functions
    @test length(random_spheres(c; radius=_ConstRadius(), n=5, rng=1)) == 5
    # :inside keeps min_gap across periodic images too
    s = random_spheres(c; radius=0.1, fraction=0.2, boundary=:inside, face_gap=0, min_gap=0.02, rng=3)
    @test min_image_gap(s, c) >= 0.02 - 1e-12
    # integer seed and Xoshiro(seed) give the same packing; edge band
    @test random_spheres(c; radius=0.08, n=10, rng=3) == random_spheres(c; radius=0.08, n=10, rng=HyperRVE.Random.Xoshiro(3))
    for x in random_spheres(c; radius=0.1, fraction=0.2, rng=5)
        t = min.(x.center, 1 .- x.center)
        @test all(abs(hypot(t[a], t[b]) - 0.1) >= 0.01 - 1e-12 for (a, b) in ((1, 2), (1, 3), (2, 3)))
    end
end
