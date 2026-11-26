# BCCZ / FCCZ presets and 2D strut lattices.

@testset "BCCZ and FCCZ" begin
    c = unitcell()
    for (k, base, nn) in ((:bccz, :bcc, 2), (:fccz, :fcc, 4))
        t = lattice_topology(k)
        @test length(t.nodes) == nn && length(t.edges) == length(lattice_topology(base).edges) + 1
        @test (1, 1, (0, 0, 1)) in t.edges
        @test relative_density(k, strut_radius(k, 0.15)) ≈ 0.15
        @test relative_density(k, 0.1) > relative_density(base, 0.1)
    end
    g = geo(RVE(c, Phase(:l, lattice(c, :bccz; density=0.12); fuse=true); matrix=nothing))
    @test volume_fraction(g, :l) ≈ 0.12 rtol = 1e-4
    g = geo(RVE(c, Phase(:l, lattice(c, :fccz; radius=0.07); fuse=true); matrix=nothing))
    @test volume_fraction(g, :l) ≈ relative_density(:fccz, 0.07) rtol = 1e-4
end

@testset "2D lattices" begin
    for (k, nn, deg, L) in ((:square, 1, 4, (1.0, 1.0)), (:triangular, 2, 6, (1.0, sqrt(3))), (:honeycomb, 4, 3, (sqrt(3), 3.0)),
                            (:kagome, 6, 4, (2.0, 2sqrt(3))), (:reentrant, 4, 3, (sqrt(3), 3.0)))
        t = lattice_topology(k)
        @test length(t.nodes) == nn && all(t.size .≈ L)
        d = zeros(Int, nn)
        for (i, j, _) in t.edges
            d[i] += 1; d[j] += 1
        end
        @test all(==(deg), d)
        lengths = [sqrt(sum(abs2, (t.nodes[j] .+ o .- t.nodes[i]) .* t.size)) for (i, j, o) in t.edges]
        @test all(l -> l ≈ 1 || (k === :reentrant && l ≈ 2), lengths)
        @test lattice_cell(k; a=0.2).size == 0.2 .* t.size
    end
    # re-entrant with angle -30° and ratio 1 is the regular honeycomb
    h, r = lattice_topology(:honeycomb), lattice_topology(:reentrant; angle=-π / 6, ratio=1.0)
    @test all(h.size .≈ r.size) && all(all(a .≈ b) for (a, b) in zip(h.nodes, r.nodes))
    # exact area of the square grid: two strips minus their overlap
    c = lattice_cell(:square; a=1.0)
    res = generate(RVE(c, Phase(:l, lattice(c, :square; thickness=0.08); fuse=true)); h=0.05)
    @test volume_fraction(res, :l) ≈ 2 * 0.08 - 0.08^2 rtol = 1e-8
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    # auxetic re-entrant honeycomb, periodic mesh, any scale and origin
    c = lattice_cell(:reentrant; a=1e-3, angle=π / 6, ratio=2.0, origin=(0.5, -1.0))
    res = generate(RVE(c, Phase(:l, lattice(c, :reentrant; thickness=8e-5); fuse=true)); h=5e-5)
    @test all(x -> x.ok, res.periodicity) && 0 < volume_fraction(res, :l) < 0.2
    for k in (:triangular, :honeycomb, :kagome)
        c = lattice_cell(k; a=0.3)
        g = geo(RVE(c, Phase(:l, lattice(c, k; thickness=0.03); fuse=true); matrix=nothing))
        @test length(g[:l].entities) >= 1
    end
    # errors
    @test_throws ArgumentError lattice(Cell((1.0, 1.0)), :honeycomb; thickness=0.05)       # wrong proportions
    @test_throws ArgumentError lattice_topology(:reentrant; angle=π / 3, ratio=0.5)        # struts would cross
    @test_throws ArgumentError lattice_topology(:hexagonal)
    @test_throws ArgumentError lattice_topology(:bcc; angle=0.1)
    @test_throws ArgumentError struts(Cell((1.0, 1.0)), [(0.0, 0.0)], [(1, 1, (1, 0))]; radius=0.05, thickness=0.1)
    @test length(struts(Cell((1.0, 1.0)), [(0.0, 0.0)], [(1, 1, (1, 0))]; thickness=0.1)) == 2
end
