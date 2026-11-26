# Strut lattices: topology of the presets, struts(), densities, meshing.

@testset "preset topologies (S1)" begin
    table = Dict(:cubic => (1, 3, 6), :bcc => (2, 8, 8), :fcc => (4, 12, nothing), :octet => (4, 24, 12),
                 :kelvin => (12, 24, 4), :diamond => (8, 16, 4))
    for (k, (nn, ne, val)) in table
        t = lattice_topology(k)
        @test length(t.nodes) == nn
        @test length(t.edges) == ne
        if val !== nothing
            deg = zeros(Int, nn)
            for (i, j, _) in t.edges
                deg[i] += 1; deg[j] += 1
            end
            @test all(==(val), deg)
        end
    end
    @test_throws ArgumentError lattice_topology(:nope)
end

@testset "struts() validation and emission (S2, S3)" begin
    c = unitcell()
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0, 0, 0)], [(1, 2)]; radius=0.05)
    @test_throws ArgumentError struts(c, [(0, 0, 0)], [(1, 1)]; radius=0.05)
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0.5, 0.5, 0.5)], [(1, 2), (2, 1)]; radius=0.05)
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0.5, 0.5, 0.5)], [(1, 3)]; radius=0.05)
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0.5, 0.5, 0.5)], [(1, 2)]; radius=0.6)
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0.5, 0.5, 0.5)], [(1, 2)]; radius=[0.1, 0.1])
    @test_throws ArgumentError lattice(c, :bcc)
    @test_throws ArgumentError lattice(Cell((1.0, 2.0, 1.0)), :bcc; density=0.1)
    # periodic cell: struts listed once (+ node spheres)
    items = lattice(c, :bcc; radius=0.05)
    @test count(x -> x isa Cylinder, items) == 8 && count(x -> x isa Sphere, items) == 2
    # non-periodic cubic cell: the complete finite truss (12 edges, 8 corners)
    np = lattice(Cell((1.0, 1.0, 1.0); periodic=false), :cubic; radius=0.05)
    @test count(x -> x isa Cylinder, np) == 12 && count(x -> x isa Sphere, np) == 8
    # node_radius options and per-strut radii
    @test isempty(filter(x -> x isa Sphere, lattice(c, :bcc; radius=0.05, node_radius=nothing)))
    @test all(x -> x.radius == 0.07, filter(x -> x isa Sphere, lattice(c, :bcc; radius=0.05, node_radius=0.07)))
    t = lattice_topology(:bcc)
    its = struts(c, t.nodes, t.edges; radius=[0.05 + 0.005k for k in 1:8])
    @test filter(x -> x isa Sphere, its)[2].radius ≈ 0.09               # :auto = largest strut at the node
    # nodes outside [0, 1) are wrapped and offsets corrected
    @test length(struts(c, [(1.0, 1.0, 1.0), (0.5, 0.5, 0.5)], [(2, 1)]; radius=0.05)) == 3
end

@testset "densities (S4, S5)" begin
    for k in (:cubic, :bcc, :fcc, :octet, :kelvin, :diamond)
        ρ = 0.15
        x = strut_radius(k, ρ)
        @test relative_density(k, x) ≈ ρ rtol = 1e-12
    end
    @test_throws ArgumentError relative_density(:octet, 0.5)
    @test_throws ArgumentError strut_radius(:octet, 0.9)
    # exact formula vs OCC volume of the fused lattice
    for (k, x) in ((:bcc, 0.08), (:octet, 0.06), (:kelvin, 0.05), (:cubic, 0.1))
        res = geo(RVE(unitcell(), Phase(:lat, lattice(unitcell(), k; radius=x); fuse=true); matrix=nothing))
        @test volume_fraction(res, :lat) ≈ relative_density(k, x) rtol = 1e-4
    end
    # density keyword and first-order upper bound
    res = geo(RVE(unitcell(), Phase(:lat, lattice(unitcell(), :bcc; density=0.2); fuse=true); matrix=nothing))
    @test volume_fraction(res, :lat) ≈ 0.2 rtol = 1e-4
    t = lattice_topology(:bcc)
    @test relative_density(unitcell(), t.nodes, t.edges; radius=0.08) > relative_density(:bcc, 0.08)
end

@testset "meshed lattices (S6, S7)" begin
    # kept cheap: one preset with the matrix, two without (the others are checked as geometry above)
    for (k, x, m) in ((:bcc, 0.08, :matrix), (:octet, 0.06, nothing), (:diamond, 0.06, nothing))
        res = generate(RVE(unitcell(), Phase(:lat, lattice(unitcell(), k; radius=x); fuse=true); matrix=m); h=0.1)
        @test all(c -> c.ok && c.npairs > 0, res.periodicity)
        @test length(res[:lat].entities) >= 1 && res.stats.nelements > 0
    end
    # finite non-periodic block, 2 x 1 x 1
    blk = generate(RVE(Cell((1.0, 1.0, 1.0); periodic=false), Phase(:lat, lattice(Cell((1.0, 1.0, 1.0); periodic=false), :bcc; radius=0.08); fuse=true);
                       matrix=nothing, tiling=(2, 1, 1)); h=0.1)
    @test length(blk[:lat].entities) == 1 && isempty(blk.periodicity)
end

@testset "offsets, faces, warnings and errors (review)" begin
    c = Cell((2.0, 1.0, 1.5); origin=(0.3, -0.2, 0.1))
    # axes are (ξj + o − ξi) .* size, also for nodes outside [0, 1)
    its = struts(c, [(1.0, 1.0, 1.0), (0.5, 0.25, 2.5)], [(2, 1, (1, 0, -1))]; radius=0.05, node_radius=nothing)
    @test only(its).axis == ((1.0, 1.0, 1.0) .+ (1, 0, -1) .- (0.5, 0.25, 2.5)) .* (2.0, 1.0, 1.5)
    # reversed strut with negated offset, and self-struts ±o, are duplicates
    @test_throws ArgumentError struts(unitcell(), [(0.1, 0.1, 0.1), (0.5, 0.5, 0.5)], [(1, 2, (1, 0, 0)), (2, 1, (-1, 0, 0))]; radius=0.05)
    @test_throws ArgumentError struts(unitcell(), [(0.1, 0.1, 0.1)], [(1, 1, (0, 0, 1)), (1, 1, (0, 0, -1))]; radius=0.05)
    @test length(struts(unitcell(), [(0.1, 0.1, 0.1), (0.5, 0.5, 0.5)], [(1, 2, (1, 0, 0)), (2, 1, (1, 0, 0))]; radius=0.05)) == 4
    # malformed input is an ArgumentError
    for e in ((1, 2, (1, 0)), (1, 2, 3), (1.5, 2), (1, 2, (0.5, 0, 0)))
        @test_throws ArgumentError struts(unitcell(), [(0.1, 0.1, 0.1), (0.5, 0.5, 0.5)], [e]; radius=0.05)
    end
    @test_throws ArgumentError struts(unitcell(), [(0.1, 0.1, 0.1), (0.5, 0.5, 0.5)], [(1, 2)]; radius=nothing)
    @test_throws ArgumentError struts(unitcell(), [(0.1, 0.1, 0.1), (0.5, 0.5, 0.5)], [(1, 2)]; radius=0.05, node_radius=:foo)
    # a strut crossing a non-periodic face appears at both faces (mixed periodicity)
    cm = Cell((1.0, 1.0, 1.0); periodic=(false, false, true))
    its = struts(cm, [(0.9, 0.0, 0.5), (0.1, 0.0, 0.5)], [(1, 2, (1, 0, 0))]; radius=0.03)
    @test count(x -> x isa Cylinder, its) == 4 && count(x -> x isa Sphere, its) == 4
    # warnings in physical units; no spurious "longer than the cell" for cubic struts
    @test_logs lattice(Cell((10.0, 10.0, 10.0)), :cubic; radius=0.6)
    @test_logs (:warn, r"shorter than their diameter") lattice(Cell((0.1, 0.1, 0.1)), :kelvin; radius=0.03)
    co = Cell((0.3, 0.3, 0.3); origin=(0.1, 0.1, 0.1))
    @test_logs RVE(co, Phase(:lat, lattice(co, :cubic; radius=0.015); fuse=true); matrix=nothing)
    # densities: exact cubic formula, fcc valences, inversion up to the range limit
    @test relative_density(:cubic, 0.1) ≈ 3π * 0.01 - 8sqrt(2) * 0.001
    deg = zeros(Int, 4)
    for (i, j, _) in lattice_topology(:fcc).edges
        deg[i] += 1; deg[j] += 1
    end
    @test sort(deg) == [4, 4, 4, 12]                    # corner to the 12 face centres
    for k in (:cubic, :bcc, :fcc, :octet, :kelvin, :diamond)
        xm = RVEs._DENSITY[k][3]
        @test relative_density(k, strut_radius(k, relative_density(k, xm))) ≈ relative_density(k, xm)
    end
    @test_throws ArgumentError lattice(unitcell(), :bcc; radius=0.05, density=0.1)
    @test_throws ArgumentError lattice(unitcell(), :bcc; density=0.9)
    # non-periodic cell: same volume fraction as the periodic lattice (complete finite truss)
    blk = Cell((1.0, 1.0, 1.0); periodic=false)
    v = volume_fraction(geo(RVE(blk, Phase(:l, lattice(blk, :cubic; radius=0.08); fuse=true); matrix=nothing)), :l)
    @test v ≈ relative_density(:cubic, 0.08) rtol = 1e-6
end

@testset "pruning of covered strut copies" begin
    slab = RVEs._cylinder_slab
    lo, hi = (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)
    @test slab((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0), 0.5, 0.05, lo, hi, 1e-9) == (0.0, 0.0)          # leaves at a corner
    t = slab((0.0, 0.0, 0.0), (1.0, -1.0, 0.0) ./ sqrt(2), 0.5, 0.05, lo, hi, 1e-9)
    @test t[1] == 0
    @test t[2] ≈ 0.05 atol = 1e-8                     # a sliver of length r
    @test all(slab((0.3, 0.3, 0.3), (0.0, 0.0, 1.0), 0.4, 0.05, lo, hi, 1e-9) .≈ (0.0, 0.4))       # inside
    @test slab((3.0, 3.0, 3.0), (0.0, 0.0, 1.0), 0.4, 0.05, lo, hi, 1e-9) === nothing                # far away
    # brute force: every sampled point of the cylinder inside the box has its axial coordinate in the slab
    rng = RVEs.Random.Xoshiro(1)
    for _ in 1:20
        a = Tuple(2 .* rand(rng, 3) .- 0.5); d = Tuple(randn(rng, 3)); d = d ./ sqrt(sum(abs2, d))
        L, r = 0.2 + rand(rng), 0.02 + 0.1rand(rng)
        sl = slab(a, d, L, r, lo, hi, 1e-9)
        for _ in 1:400
            tt = L * rand(rng); w = Tuple(randn(rng, 3)); w = w .- sum(w .* d) .* d; w = w ./ sqrt(sum(abs2, w)) .* (r * sqrt(rand(rng)))
            q = a .+ tt .* d .+ w
            if all(0 .<= q .<= 1)
                @test sl !== nothing && sl[1] - 1e-6 <= tt <= sl[2] + 1e-6
            end
        end
    end
    # nearly axial struts close to a face (the slice distance must not report a miss)
    sd = RVEs._slice_dist
    ε = 1e-3
    @test sd((0.5, 1 + 2e-5, 0.5), (ε, sqrt(1 - ε^2), 0.0), lo, hi) < 0.03
    for _ in 1:200
        e = 10.0^(-2 - 3rand(rng)); dd = (e * randn(rng), 1.0, e * randn(rng)); dd = dd ./ sqrt(sum(abs2, dd))
        a = (rand(rng), 1 + 10.0^(-5 + 2rand(rng)), rand(rng)) .- 0.3 .* dd
        L, r = 0.6, 0.05
        sl = slab(a, dd, L, r, lo, hi, 1e-9)
        for _ in 1:200
            tt = L * rand(rng); w = Tuple(randn(rng, 3)); w = w .- sum(w .* dd) .* dd; w = w ./ sqrt(sum(abs2, w)) .* (r * sqrt(rand(rng)))
            q = a .+ tt .* dd .+ w
            if all(0 .<= q .<= 1)
                @test sl !== nothing && sl[1] - 1e-9 <= tt <= sl[2] + 1e-9
            end
        end
    end
    cyl = Cylinder((0.5, 1e-3, 0.5), 0.5 .* (8e-3, sqrt(1 - 6.4e-5), 0.0), 0.2)
    g = geo(RVE(unitcell(), Phase(:l, [cyl]; fuse=true); matrix=nothing))
    @test volume_fraction(g, :l) > 0
    # exactness: the same body with and without pruning
    c = unitcell()
    for (k, x) in ((:bcc, 0.08), (:cubic, 0.1))
        v = map((true, false)) do p
            RVEs._PRUNE[] = p
            try
                volume_fraction(geo(RVE(c, Phase(:l, lattice(c, k; radius=x); fuse=true); matrix=nothing)), :l)
            finally
                RVEs._PRUNE[] = true
            end
        end
        @test v[1] ≈ v[2] rtol = 1e-7
    end
    # a free direction and no node spheres
    cf = Cell((1.0, 1.0, 1.0); periodic=(true, true, false))
    @test volume_fraction(geo(RVE(cf, Phase(:l, lattice(cf, :bcc; radius=0.08); fuse=true); matrix=nothing)), :l) ≈
          relative_density(:bcc, 0.08) rtol = 1e-3
    @test volume_fraction(geo(RVE(c, Phase(:l, lattice(c, :cubic; radius=0.1, node_radius=nothing); fuse=true); matrix=nothing)), :l) > 0
end
