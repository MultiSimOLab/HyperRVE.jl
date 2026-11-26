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
    # 27 pieces: the copy at the interior lattice point (1, 1, 1) of the supercell is one whole piece;
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

@testset "more supercell cases" begin
    # mixed periodicity, meshed, independently checked
    cell = Cell((1.0, 1.0, 1.0); origin=(0.5, -0.5, 2.0), periodic=(true, false, true))
    rve = RVE(cell, Phase(:p, [Sphere(cellpoint(cell, (0, 0.5, 0)), 0.2), Fibre(cellpoint(cell, (0.5, 0.5, 0.5)), (0, 1, 0), 0.1)]);
              tiling=(2, 2, 1))
    generate(rve; h=0.2, output=joinpath(OUTDIR, "super_mixed.msh")) do res
        @test [c.direction for c in res.periodicity] == [1, 3]
        @test all(c -> c.ok && c.npairs > 0, res.periodicity)
        @test brute_periodic_match("xmin", "xmax", 1, 2.0)
        @test brute_periodic_match("zmin", "zmax", 3, 1.0)
        @test mesh_volume() ≈ 4 rtol = 1e-2
    end
    # fibre volume scales with the tiling
    fib = RVE(unitcell(), Phase(:f, Fibre((0.5, 0.0, 0.5), (0, 1, 0), 0.1)); tiling=(2, 3, 1))
    @test volume(geo(fib), :f) ≈ 2 * 3 * π * 0.01 rtol = 1e-6
    # item larger than the cell
    big = RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.6)))
    @test volume_fraction(geo(supercell(big, 2)), :p) ≈ volume_fraction(geo(big), :p) rtol = 1e-6
    # no matrix, with refinement
    lat = RVE(unitcell(), Phase(:s, Sphere((0.0, 0.0, 0.0), 0.3); refine=InterfaceRefinement(level=2)); matrix=nothing)
    res = generate(supercell(lat, 2); h=0.2)
    @test all(c -> c.ok, res.periodicity)
    @test phase_names(res) == ["s"]
    @test volume(res, :s) ≈ 8 * sphere_vol(0.3) rtol = 1e-6
    # reopened file: physical coordinates of the supercell
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(joinpath(OUTDIR, "super_mixed.msh"))
        x = reshape(gmsh.model.mesh.getNodes()[2], 3, :)
        @test minimum(x[1, :]) ≈ 0.5 && maximum(x[1, :]) ≈ 2.5
        @test minimum(x[2, :]) ≈ -0.5 && maximum(x[2, :]) ≈ 1.5
    finally
        gmsh.finalize()
    end
    # a rotated ellipsoid cut by a face: OCC volumes of B-spline pieces are ~1e-2 accurate
    e = Rotate(Ellipsoid((0.0, 0.5, 0.5), (0.3, 0.12, 0.12)); axis=(1, 1, 0), angle=π / 4)
    u = RVE(unitcell(), Phase(:e, e))
    @test volume(geo(supercell(u, (2, 1, 1))), :e) ≈ 2 * 4π * 0.3 * 0.12 * 0.12 / 3 rtol = 2e-2
end

@testset "non-periodic tiling places items per copy (documented)" begin
    cell = Cell((1.0, 1.0, 1.0); periodic=(false, true, true))
    outside = RVE(cell, Phase(:p, Sphere((-0.5, 0.5, 0.5), 0.2)); tiling=(2, 1, 1))
    @test volume(geo(outside), :p) ≈ sphere_vol(0.2) rtol = 1e-6
end
