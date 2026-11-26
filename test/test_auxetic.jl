# Chiral lattices and rotating units (auxetic 2D cells).

@testset "chiral lattices" begin
    for (kind, r, t) in ((:hexachiral, 0.2, 0.04), (:tetrachiral, 0.25, 0.05), (:trichiral, 0.2, 0.05),
                         (:antitetrachiral, 0.25, 0.05), (:antitrichiral, 0.2, 0.05))
        c = lattice_cell(kind; a=1.0)
        ρ = relative_density(kind; ring_radius=r, thickness=t)
        for mirror in (false, true)
            g = geo(RVE(c, Phase(:l, chiral(c, kind; ring_radius=r, thickness=t, mirror); fuse=true); matrix=nothing))
            @test volume_fraction(g, :l) ≈ ρ rtol = 1e-9
        end
        g = geo(RVE(c, Phase(:l, chiral(c, kind; ring_radius=r, thickness=t, filled=true); fuse=true); matrix=nothing))
        @test volume_fraction(g, :l) ≈ relative_density(kind; ring_radius=r, thickness=t, filled=true) rtol = 1e-9
        # small cell away from the origin: periodic mesh in one piece
        cs = lattice_cell(kind; a=1e-3, origin=(0.5, -1.0))
        nc = Ref(0)
        res = generate(RVE(cs, Phase(:l, chiral(cs, kind; ring_radius=r * 1e-3, thickness=t * 1e-3); fuse=true));
                       h=0.05e-3) do _
            nc[] = periodic_components_2d(cs, "l")
        end
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
        @test volume_fraction(res, :l) ≈ ρ rtol = 1e-6
        @test nc[] == 1
        @test res.stats.min_quality > 0
    end
    c = lattice_cell(:hexachiral)
    res = generate(RVE(c, Phase(:l, chiral(c, :hexachiral; ring_radius=0.2, thickness=0.04); fuse=true)); h=0.05) do _
        @test brute_periodic_match_2d("xmin", "xmax", 1, c.size[1])
        @test brute_periodic_match_2d("ymin", "ymax", 2, c.size[2])
    end
    # handedness: the mirror image is a different geometry with the same area
    a1 = chiral(c, :hexachiral; ring_radius=0.2, thickness=0.04)
    a2 = chiral(c, :hexachiral; ring_radius=0.2, thickness=0.04, mirror=true)
    @test a1[end].vertices != a2[end].vertices
    @test lattice_cell(:tetrachiral; a=0.5).size == (0.5, 0.5)
    @test all(lattice_cell(:hexachiral; a=2.0).size .≈ (2.0, 2sqrt(3)))
end

@testset "rotating units" begin
    for θ in (π / 12, π / 6, π / 4, π / 3)
        c = lattice_cell(:rotating_squares; a=1.0, angle=θ)
        @test all(c.size .≈ 2 * (cos(θ) + sin(θ)))
        g = geo(RVE(c, Phase(:l, rotating_squares(c; angle=θ, hinge=0.1); fuse=true); matrix=nothing))
        @test volume_fraction(g, :l) ≈ relative_density(:rotating_squares; angle=θ, hinge=0.1) rtol = 1e-9
    end
    @test relative_density(:rotating_squares; angle=π / 8, hinge=0.1) ≈ relative_density(:rotating_squares; angle=3π / 8, hinge=0.1)
    for θ in (π / 6, π / 3)
        c = lattice_cell(:rotating_triangles; a=1.0, angle=θ)
        g = geo(RVE(c, Phase(:l, rotating_triangles(c; angle=θ, hinge=0.1); fuse=true); matrix=nothing))
        @test volume_fraction(g, :l) ≈ relative_density(:rotating_triangles; angle=θ, hinge=0.1) rtol = 1e-9
    end
    # control for the component counter: squares that only touch at corners are separate pieces
    θ = π / 6
    c = lattice_cell(:rotating_squares; angle=θ)
    sq = [it for it in rotating_squares(c; angle=θ, hinge=0.1) if it isa Polygon]
    nc = Ref(0)
    generate(RVE(c, Phase(:l, sq; fuse=true)); h=0.1) do _
        nc[] = periodic_components_2d(c, "l")
    end
    @test nc[] == 4
    for (f, kind) in ((rotating_squares, :rotating_squares), (rotating_triangles, :rotating_triangles))
        θ = deg2rad(15)
        c = lattice_cell(kind; a=1e-3, angle=θ, origin=(0.5, -1.0))
        nc = Ref(0)
        res = generate(RVE(c, Phase(:l, f(c; angle=θ, hinge=0.05e-3); fuse=true)); h=0.07e-3) do _
            nc[] = periodic_components_2d(c, "l")
        end
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
        @test volume_fraction(res, :l) ≈ relative_density(kind; angle=θ, hinge=0.05) rtol = 1e-6
        @test nc[] == 1
        @test res.stats.min_quality > 0
    end
end

@testset "auxetic validation" begin
    c = lattice_cell(:hexachiral)
    @test_throws ArgumentError chiral(c, :pentachiral; ring_radius=0.2, thickness=0.04)
    @test_throws ArgumentError chiral(Cell((1.0, 1.0)), :hexachiral; ring_radius=0.2, thickness=0.04)       # proportions
    @test_throws ArgumentError chiral(Cell(c.size; periodic=(true, false)), :hexachiral; ring_radius=0.2, thickness=0.04)
    @test_throws ArgumentError chiral(c, :hexachiral; ring_radius=0.2)
    @test_throws ArgumentError chiral(c, :hexachiral; ring_radius=0.2, thickness=0.5)                      # t >= 2r
    @test_throws ArgumentError chiral(c, :hexachiral; ring_radius=0.45, thickness=0.2)                     # rings touch
    @test_throws ArgumentError relative_density(:hexachiral; ring_radius=0.2, thickness=0.15)              # beyond exact range
    @test_throws ArgumentError lattice(c, :hexachiral; thickness=0.05)
    @test_throws ArgumentError lattice_cell(:hexachiral; angle=0.1)
    cs = lattice_cell(:rotating_squares; angle=π / 6)
    @test_throws ArgumentError rotating_squares(cs; angle=π / 6)
    @test_throws ArgumentError rotating_squares(cs; angle=π / 2, hinge=0.1)
    @test_throws ArgumentError rotating_squares(cs; angle=π / 6, hinge=0.0)
    @test_throws ArgumentError rotating_squares(cs; angle=π / 6, hinge=1.1)                                # wider than 2a sin θ
    @test_throws ArgumentError rotating_squares(lattice_cell(:hexachiral); angle=π / 6, hinge=0.1)
    @test_throws ArgumentError lattice_cell(:rotating_squares)
    @test_throws ArgumentError lattice(cs, :rotating_squares; thickness=0.05)
    ct = lattice_cell(:rotating_triangles; angle=π / 6)
    @test_throws ArgumentError rotating_triangles(ct; angle=π / 6, hinge=1.2)
    @test_throws ArgumentError rotating_triangles(ct; angle=2π / 3, hinge=0.1)
end

@testset "auxetic review fixes" begin
    # filled: exact up to the overlap of neighbouring ligaments, error beyond
    for (kind, r, tmax) in ((:tetrachiral, 0.1, 0.2), (:antitetrachiral, 0.1, 0.2), (:trichiral, 0.05, 0.3), (:antitrichiral, 0.05, 0.3))
        c = lattice_cell(kind)
        t = 0.97 * tmax
        g = geo(RVE(c, Phase(:l, chiral(c, kind; ring_radius=r, thickness=t, filled=true); fuse=true); matrix=nothing))
        @test volume_fraction(g, :l) ≈ relative_density(kind; ring_radius=r, thickness=t, filled=true) rtol = 1e-9
        @test_throws ArgumentError relative_density(kind; ring_radius=r, thickness=1.1tmax, filled=true)
    end
    c = lattice_cell(:hexachiral)
    g = geo(RVE(c, Phase(:l, chiral(c, :hexachiral; ring_radius=0.2, thickness=0.133); fuse=true); matrix=nothing))
    @test volume_fraction(g, :l) ≈ relative_density(:hexachiral; ring_radius=0.2, thickness=0.133) rtol = 1e-9
    # anti-chiral with filled rings and t = 2r: ligament edges on the faces
    for kind in (:antitetrachiral, :antitrichiral)
        @test_throws ArgumentError chiral(lattice_cell(kind), kind; ring_radius=0.1, thickness=0.2, filled=true)
        c = lattice_cell(kind)
        res = generate(RVE(c, Phase(:l, chiral(c, kind; ring_radius=0.1, thickness=0.21, filled=true); fuse=true)); h=0.08)
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    end
    # same lengths as the geometry with a
    @test relative_density(:hexachiral; ring_radius=2e-4, thickness=4e-5, a=1e-3) ≈ relative_density(:hexachiral; ring_radius=0.2, thickness=0.04)
    @test relative_density(:rotating_squares; angle=0.3, hinge=1e-4, a=1e-3) ≈ relative_density(:rotating_squares; angle=0.3, hinge=0.1)
    # foreign keywords and wrong calls
    @test_throws ArgumentError relative_density(:rotating_squares; angle=0.3, hinge=0.1, thickness=0.1)
    @test_throws ArgumentError relative_density(:hexachiral; ring_radius=0.2, thickness=0.04, angle=0.3)
    @test_throws ArgumentError relative_density(:rotating_squares; angle=0.3)
    @test_throws ArgumentError relative_density(:rotating_triangles; angle=0.3, hinge=2.0)
    @test_throws ArgumentError relative_density(:square; ring_radius=0.2, thickness=0.04)
    @test_throws ArgumentError relative_density(:hexachiral, 0.1)
    @test_throws ArgumentError lattice_topology(:hexachiral)
    @test_throws ArgumentError lattice_cell(:pentachiral)
    @test_throws ArgumentError lattice_cell(:rotating_squares; angle="30")
    @test_throws ArgumentError lattice_cell(:rotating_squares; angle=0.0)
    ct = lattice_cell(:rotating_triangles; angle=0.5)
    @test_throws ArgumentError rotating_triangles(Cell(ct.size; periodic=(true, false)); angle=0.5, hinge=0.1)
    @test_throws ArgumentError rotating_triangles(Cell((1.0, 1.0)); angle=0.5, hinge=0.1)
    # rotating units near the limits, and triangles beyond the kagome angle
    for (kind, f, θ, w) in ((:rotating_squares, rotating_squares, 0.2, 0.999 * 2sin(0.2)),
                            (:rotating_triangles, rotating_triangles, 1.6, 0.999 * 2sin(2π / 3 - 1.6)),
                            (:rotating_triangles, rotating_triangles, 1.2, 0.1))
        c = lattice_cell(kind; angle=θ)
        g = geo(RVE(c, Phase(:l, f(c; angle=θ, hinge=w); fuse=true); matrix=nothing))
        @test volume_fraction(g, :l) ≈ relative_density(kind; angle=θ, hinge=w) rtol = 1e-9
    end
    # cells of every kind scale with a
    for kind in (:hexachiral, :tetrachiral, :trichiral, :antitetrachiral, :antitrichiral)
        @test all(lattice_cell(kind; a=0.3).size .≈ 0.3 .* lattice_cell(kind).size)
    end
    for kind in (:rotating_squares, :rotating_triangles)
        @test all(lattice_cell(kind; a=0.3, angle=0.4).size .≈ 0.3 .* lattice_cell(kind; angle=0.4).size)
    end
    # handedness: the mirror of a chiral lattice is its reflection, the anti-chiral ones are achiral
    c = lattice_cell(:tetrachiral)
    cy = cellpoint(c, lattice_topology(:square).nodes[1])[2]
    wrap(p) = round.(mod.(p, 1.0) .% 1.0; digits=8)
    pts(items) = sort(unique([wrap(v) for it in items if it isa Polygon for v in it.vertices]))
    a = chiral(c, :tetrachiral; ring_radius=0.25, thickness=0.05)
    b = chiral(c, :tetrachiral; ring_radius=0.25, thickness=0.05, mirror=true)
    refl = [Polygon([(v[1], 2cy - v[2]) for v in it.vertices]) for it in a if it isa Polygon]
    @test pts(refl) == pts(b) && pts(a) != pts(b)
    # periodic matching on an anti-chiral and a rotating cell
    for (c, items) in ((lattice_cell(:antitetrachiral), chiral(lattice_cell(:antitetrachiral), :antitetrachiral; ring_radius=0.25, thickness=0.05)),
                       (lattice_cell(:rotating_squares; angle=0.4), rotating_squares(lattice_cell(:rotating_squares; angle=0.4); angle=0.4, hinge=0.1)))
        generate(RVE(c, Phase(:l, items; fuse=true)); h=0.08) do _
            @test brute_periodic_match_2d("xmin", "xmax", 1, c.size[1])
            @test brute_periodic_match_2d("ymin", "ymax", 2, c.size[2])
        end
    end
    # control for the counter with triangles that only touch at corners
    θ = 0.5
    c = lattice_cell(:rotating_triangles; angle=θ)
    tri = [it for it in rotating_triangles(c; angle=θ, hinge=0.1) if it isa Polygon]
    nc = Ref(0)
    generate(RVE(c, Phase(:l, tri; fuse=true)); h=0.1) do _
        nc[] = periodic_components_2d(c, "l")
    end
    @test nc[] == 4
    # finite samples: supercell of an auxetic cell
    c = lattice_cell(:hexachiral)
    rve = supercell(RVE(c, Phase(:l, chiral(c, :hexachiral; ring_radius=0.2, thickness=0.04); fuse=true)), (2, 1))
    res = generate(rve; h=0.06)
    @test volume_fraction(res, :l) ≈ relative_density(:hexachiral; ring_radius=0.2, thickness=0.04) rtol = 1e-6
end

