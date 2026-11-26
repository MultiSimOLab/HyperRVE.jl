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
