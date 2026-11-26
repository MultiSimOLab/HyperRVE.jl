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
