# Periodic Voronoi tessellations and polycrystals.

"V − E + F − C of a fully periodic 3D tessellation (0 on the 3-torus)."
function voronoi_euler(t)
    E = Set{Any}(); F = Set{Any}()
    for a in eachindex(t.loops), (f, (b, k)) in zip(t.loops[a], t.neighbours[a])
        push!(F, min((a, b, k), (b, a, .-k)))
        for i in eachindex(f)
            (v1, o1), (v2, o2) = f[i], f[mod1(i + 1, length(f))]
            push!(E, min((v1, v2, o2 .- o1), (v2, v1, o1 .- o2)))
        end
    end
    return length(t.vertices) - length(E) + length(F) - length(t.loops)
end

@testset "tessellation invariants" begin
    c = unitcell()
    for (n, s) in ((1, 1), (2, 1), (8, 2), (30, 3), (60, 4))
        t = voronoi(c; n, rng=s)
        @test sum(grain_volumes(t)) ≈ 1 atol = 1e-12
        n >= 2 && @test voronoi_euler(t) == 0
    end
    t = voronoi(c; n=30, min_distance=0, rng=5)                       # Poisson
    @test sum(grain_volumes(t)) ≈ 1 atol = 1e-12
    @test voronoi_euler(t) == 0
    # degenerate seeds: cubes, truncated octahedra, rhombic dodecahedra
    nf(t) = sort(unique(length.(t.loops)))
    @test nf(voronoi(c; seeds=[(x, y, z) for x in (0.25, 0.75) for y in (0.25, 0.75) for z in (0.25, 0.75)])) == [6]
    @test nf(voronoi(c; seeds=[(0, 0, 0), (0.5, 0.5, 0.5)])) == [14]
    t = voronoi(c; seeds=[(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)])
    @test nf(t) == [12] && all(grain_volumes(t) .≈ 0.25)
    # non-cubic cell with origin, mixed periodicity, 2D
    cb = Cell((2.0, 1.0, 0.5); origin=(1.0, -1.0, 3.0), periodic=(true, false, true))
    t = voronoi(cb; n=12, rng=7)
    @test sum(grain_volumes(t)) ≈ 1.0 atol = 1e-12
    @test all(g -> all(v -> -1e-12 <= v[2] + 1.0 <= 1.0 + 1e-12, g.vertices), grains(t))   # no wrap along y
    t2 = voronoi(Cell((1.0, 2.0)); n=20, rng=8)
    @test sum(grain_volumes(t2)) ≈ 2 atol = 1e-12
    @test grain_volumes(voronoi(Cell((1.0, 1.0)); seeds=[(0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75)])) ≈ fill(0.25, 4)
    # seeds, Lloyd, reproducibility, errors
    @test voronoi(c; n=10, rng=3).seeds == voronoi(c; n=10, rng=3).seeds
    spread(t) = (v = grain_volumes(t); maximum(v) / minimum(v))
    @test spread(voronoi(c; n=20, lloyd=5, rng=9)) < spread(voronoi(c; n=20, rng=9))
    @test all(voronoi(c; seeds=[(1.25, -0.5, 0.5), (0.5, 0.5, 0.75)]).seeds[1] .≈ (0.25, 0.5, 0.5))
    @test_throws ArgumentError voronoi(c)
    @test_throws ArgumentError voronoi(c; n=3, seeds=[(0.1, 0.1, 0.1)])
    @test_throws ArgumentError voronoi(c; seeds=[(0.1, 0.1, 0.1), (1.1, 0.1, 0.1)])                    # coincide
    @test_throws ArgumentError voronoi(Cell((1.0, 1.0, 1.0); periodic=false); seeds=[(0.1, 0.1, 1.5)])
    @test_throws ArgumentError voronoi(c; n=1000, min_distance=0.5)
    @test_throws ArgumentError voronoi(c; n=0)
    @test_throws ArgumentError grain_phases(voronoi(c; n=4, rng=1); groups=[[1, 2], [2, 3, 4]])
end

@testset "polycrystals" begin
    c = unitcell()
    t = voronoi(c; n=8, rng=8)
    rve = RVE(c, grain_phases(t); matrix=nothing)
    g = geo(rve)
    @test [volume(g, Symbol("grain_", i)) for i in 1:8] ≈ grain_volumes(t) atol = 1e-9
    res = generate(rve; h=0.1)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test length(interface_names(res)) > 8 && all(startswith("interface_grain_"), interface_names(res))
    # grain groups as fused phases, free faces
    t = voronoi(Cell((1.0, 1.0, 0.5); periodic=(true, true, false)); n=10, rng=3)
    res = generate(RVE(t.cell, grain_phases(t; groups=[mod1(i, 2) for i in 1:10]); matrix=nothing); h=0.1)
    @test sum(p.fraction for p in res.phases) ≈ 1 atol = 1e-8
    @test all(x -> x.ok, res.periodicity)
    # 2D
    c2 = Cell((1.0, 1.0))
    t = voronoi(c2; n=12, rng=1)
    res = generate(RVE(c2, grain_phases(t); matrix=nothing); h=0.04)
    @test [volume(res, Symbol("grain_", i)) for i in 1:12] ≈ grain_volumes(t) atol = 1e-9
    @test all(x -> x.ok, res.periodicity)
end

@testset "review fixes" begin
    c = unitcell()
    # perturbed regular lattices: the weld merges near-coincident vertices consistently
    lattices = Dict(:cubic2 => [(x, y, z) for x in (0.25, 0.75) for y in (0.25, 0.75) for z in (0.25, 0.75)],
                    :cubic3 => [(x, y, z) for x in (1, 3, 5) ./ 6 for y in (1, 3, 5) ./ 6 for z in (1, 3, 5) ./ 6],
                    :fcc => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)])
    rng = RVEs.Random.Xoshiro(1)
    nfail = 0
    for (_, S) in lattices, a in (1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-3), _ in 1:3
        seeds = [p .+ a .* (2 .* (rand(rng), rand(rng), rand(rng)) .- 1) for p in S]
        try
            t = voronoi(c; seeds)
            g = grains(t)
            abs(sum(grain_volumes(t)) - 1) < 1e-9 || (nfail += 1)
        catch e
            nfail += 1
        end
    end
    @test nfail == 0
    for a in (1e-8, 1e-7, 1e-6, 1e-4), _ in 1:3
        S = [(x, y) .+ a .* (2 .* (rand(rng), rand(rng)) .- 1) for x in (1, 3, 5, 7) ./ 8 for y in (1, 3, 5, 7) ./ 8]
        t = voronoi(Cell((1.0, 1.0)); seeds=S)
        @test sum(grain_volumes(t)) ≈ 1 atol = 1e-9
        @test length(grains(t)) == 16
    end
    # a periodic tessellation cut by free faces fills the cell
    cn = Cell((1.0, 1.0, 1.0); periodic=false)
    t = voronoi(cn; n=6, rng=1, periodic=true)
    g = geo(RVE(cn, grain_phases(t); matrix=nothing))
    @test sum(p.fraction for p in g.phases) ≈ 1 atol = 1e-8
    @test [volume(g, Symbol("grain_", i)) for i in 1:6] ≈ grain_volumes(t) atol = 1e-8
    c2 = Cell((1.0, 1.0); periodic=(true, false))
    t = voronoi(c2; n=8, rng=2, periodic=true)
    res = generate(RVE(c2, grain_phases(t); matrix=nothing); h=0.05)
    @test sum(p.fraction for p in res.phases) ≈ 1 atol = 1e-8
    # :auto distance in slender cells and with one seed
    @test length(voronoi(Cell((2.0, 1.0, 0.5)); n=1).seeds) == 1
    @test length(voronoi(Cell((1.0, 1.0, 0.05)); n=3, rng=1).seeds) == 3
    @test length(voronoi(Cell((10.0, 1.0)); n=2, rng=1).seeds) == 2
    # group labels
    t = voronoi(c; n=4, rng=1)
    @test [p.name for p in grain_phases(t; groups=[1, 3, 1, 3])] == ["grain_1", "grain_3"]
    @test [p.name for p in grain_phases(t; groups=[:a, :b, :a, :b])] == ["a", "b"]
    @test [p.name for p in grain_phases(t; groups=Any[2, 1, 2, 1], prefix="ph")] == ["ph_1", "ph_2"]
    @test all(p -> p.fuse, grain_phases(t; groups=[1, 2, 1, 2]))
    @test grain_phases(t; refine=0.01)[1].refine !== nothing
    @test_throws ArgumentError grain_phases(t; groups=[1, 2, 1])
    @test_throws ArgumentError grain_phases(t; groups=[true, false, true, false])
    # input validation
    @test_throws ArgumentError voronoi(c; seeds=[(0.5, 0.5, 0.5), (0.5, 0.5, 0.5 + 3e-8)])
    @test_throws ArgumentError voronoi(c; n=4, periodic=(true, false))
    @test_throws ArgumentError voronoi(c; n=true)
end
