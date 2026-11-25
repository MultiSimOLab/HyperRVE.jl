# Supercells: the phases of one cell repeated n1 x n2 x n3 times.

@testset "volumes scale with the number of cells (U1)" begin
    unit = RVE(unitcell(), Phase(:p, [Sphere((0.5, 0.5, 0.5), 0.2), Sphere((0.0, 0.0, 0.0), 0.15)]))
    v1 = geo(unit)
    for n in ((2, 2, 2), (3, 1, 2))
        vn = geo(supercell(unit, n))
        @test volume(vn) ≈ prod(n) rtol = 1e-12
        @test volume(vn, :p) ≈ prod(n) * volume(v1, :p) rtol = 1e-6
        @test volume_fraction(vn, :p) ≈ volume_fraction(v1, :p) rtol = 1e-6
    end
end

@testset "corner sphere in a periodic 2x2x2 supercell (U2)" begin
    r = 0.2
    res = geo(supercell(RVE(unitcell(), Phase(:p, Sphere((0, 0, 0), r))), 2))
    @test volume(res, :p) ≈ 8 * sphere_vol(r) rtol = 1e-6
    # the copy at the interior lattice point (1, 1, 1) of the supercell is one whole piece;
    # the others are split by the supercell faces into 2, 4 or 8 pieces
    @test length(res[:p].entities) == 1 + 3 * 2 + 3 * 4 + 8
end

@testset "non-periodic direction (U3)" begin
    cell = Cell((1.0, 1.0, 1.0); periodic=(false, true, true))
    # crosses the outer xmin face (cut) and, in copy 2, the internal cell boundary (whole)
    res = geo(RVE(cell, Phase(:p, Sphere((0.0, 0.5, 0.5), 0.2)); tiling=(2, 1, 1)))
    @test volume(res, :p) ≈ 1.5 * sphere_vol(0.2) rtol = 1e-6
    res2 = geo(RVE(cell, Phase(:p, Sphere((1.0, 0.5, 0.5), 0.2)); tiling=(2, 1, 1)))
    @test volume(res2, :p) ≈ 1.5 * sphere_vol(0.2) rtol = 1e-6
end

@testset "meshed supercells are periodic (U4)" begin
    e = Rotate(Ellipsoid((0.5, 0.5, 0.5), (0.3, 0.12, 0.12)); axis=(1, 1, 0), angle=π / 4)
    unit = RVE(Cell((1.0, 1.0, 1.0); origin=(2.0, 0.0, -1.0)), Phase(:e, e))
    g = geo(supercell(unit, (3, 3, 3)))
    @test volume(g, :e) ≈ 27 * volume(geo(unit), :e) rtol = 1e-6
    res = generate(supercell(RVE(unitcell(), Phase(:p, [Sphere((0, 0, 0), 0.25), Fibre((0.5, 0.5, 0.5), (0, 0, 1), 0.1)])),
                             (2, 1, 1)); h=0.15)
    @test all(c -> c.ok && c.npairs > 0, res.periodicity)
    @test volume(res) ≈ 2
    names = group_names(res; dim=2)
    @test "xmin" in names && "xmax" in names
    # the anisotropic domain: cells of 2 x 1 x 0.5 tiled (1, 2, 4) -> 2 x 2 x 2 domain
    c2 = Cell((2.0, 1.0, 0.5); origin=(-1.0, 3.0, 0.0))
    r2 = generate(RVE(c2, Phase(:p, Sphere(cellpoint(c2, (0, 0, 0)), 0.2)); tiling=(1, 2, 4)); h=0.25)
    @test volume(r2) ≈ 8
    @test volume(r2, :p) ≈ 8 * sphere_vol(0.2) rtol = 1e-6
    @test all(c -> c.ok && c.npairs > 0, r2.periodicity)
end
