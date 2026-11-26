# Geometry invariants (geometry_only = true: no mesh, exact OCC volumes).

r = 0.2
@testset "centred sphere (G1)" begin
    for per in (true, false)
        res = geo(RVE(unitcell(periodic=per), Phase(:p, Sphere((0.5, 0.5, 0.5), r))))
        @test volume(res, :p) ≈ sphere_vol(r) rtol = 1e-6
        @test sum(last, volume_fractions(res)) ≈ 1 rtol = 1e-9
        @test phase_names(res) == ["matrix", "p"]
        @test res.stats === nothing && isempty(res.periodicity)
    end
end

@testset "sphere on a face, edge and corner (G2)" begin
    for c in ((0.0, 0.5, 0.5), (0.0, 0.0, 0.5), (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
        res = geo(RVE(unitcell(), Phase(:p, Sphere(c, r))))
        @test volume(res, :p) ≈ sphere_vol(r) rtol = 1e-6
    end
    # entities: a corner sphere is split in 8 pieces
    @test length(geo(RVE(unitcell(), Phase(:p, Sphere((0, 0, 0), r))))[:p].entities) == 8
end

@testset "sphere larger than the cell (G3)" begin
    R = 0.6
    a = geo(RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), R))))
    b = geo(RVE(unitcell(), Phase(:p, Sphere((0.0, 0.0, 0.0), R))))
    @test 0.5 < volume_fraction(a, :p) < sphere_vol(R)
    @test volume(a, :p) ≈ volume(b, :p) rtol = 1e-6
end

@testset "BCC and FCC (G4)" begin
    bcc = geo(RVE(unitcell(), Phase(:p, [Sphere((0.5, 0.5, 0.5), 0.25), Sphere((1, 1, 1), 0.25)])))
    @test volume(bcc, :p) ≈ 2 * sphere_vol(0.25) rtol = 1e-6
    fcc = geo(RVE(unitcell(), Phase(:p, [Sphere((0, 0.5, 0.5), 0.25), Sphere((0.5, 0, 0.5), 0.25),
                                       Sphere((0.5, 0.5, 0), 0.25), Sphere((1, 1, 1), 0.25)])))
    @test volume(fcc, :p) ≈ 4 * sphere_vol(0.25) rtol = 1e-6
end

@testset "fibres (G5)" begin
    rf = 0.1
    centre = geo(RVE(unitcell(), Phase(:f, Fibre((0.5, 0.5, 0.5), (1, 0, 0), rf))))
    @test volume(centre, :f) ≈ π * rf^2 rtol = 1e-6
    edge = geo(RVE(unitcell(), Phase(:f, Fibre((0.3, 0.0, 1.0), (1, 0, 0), rf))))
    @test volume(edge, :f) ≈ π * rf^2 rtol = 1e-6
    # 4x4 array with the last row/column on the opposite faces: 3x3 fibres per cell
    st = 1 / 3
    rr = 0.9 * st / 2
    arr = geo(RVE(unitcell(), Phase(:f, [Fibre((0.5, i * st, j * st), (1, 0, 0), rr) for i in 0:3 for j in 0:3])))
    @test volume(arr, :f) ≈ 9π * rr^2 rtol = 1e-6
    # three orthogonal fibres as three items of one phase
    c = (0.5, 0.5, 0.5)
    three = geo(RVE(unitcell(), Phase(:f, [Fibre(c, (1, 0, 0), rf), Fibre(c, (0, 1, 0), rf), Fibre(c, (0, 0, 1), rf)])))
    @test π * rf^2 < volume(three, :f) < 3π * rf^2
    # a non-unit anisotropic cell: fibre along y in a 2 x 0.5 x 1 cell
    cell = Cell((2.0, 0.5, 1.0); origin=(-1.0, 3.0, 0.0))
    ff = geo(RVE(cell, Phase(:f, Fibre(cellpoint(cell, (0, 0, 1)), (0, 1, 0), rf))))
    @test volume(ff, :f) ≈ π * rf^2 * 0.5 rtol = 1e-6
end

@testset "rotated ellipsoid (G6)" begin
    e = Rotate(Rotate(Ellipsoid((0.5, 0.5, 0.5), (0.25, 0.1, 0.1)); axis=(0, 1, 0), angle=π / 4);
               axis=(0, 0, 1), angle=π / 4)
    res = geo(RVE(unitcell(), Phase(:e, e)))
    # the ellipsoid is a dilated sphere (B-spline in OCC): volume accurate to ~1e-4
    @test volume(res, :e) ≈ 4π * 0.25 * 0.1 * 0.1 / 3 rtol = 1e-3
    # crossing the corner: still complete
    res2 = geo(RVE(unitcell(), Phase(:e, Rotate(Ellipsoid((1, 1, 1), (0.25, 0.1, 0.1)); angle=0.4, axis=(1, 1, 0)))))
    @test volume(res2, :e) ≈ 4π * 0.25 * 0.1 * 0.1 / 3 rtol = 1e-3
end

@testset "booleans and transforms (G7)" begin
    c = (0.5, 0.5, 0.5)
    L = 0.38
    cross = Fuse(Fibre(c, (1, 0, 0), 0.1), Fibre(c, (0, 1, 0), 0.1), Fibre(c, (0, 0, 1), 0.1))
    core = Intersect(Sphere(c, 0.25), Box(center=c, size=(L, L, L)))
    v = volume(geo(RVE(unitcell(), Phase(:i, Cut(core, cross)))), :i)
    vcore = volume(geo(RVE(unitcell(), Phase(:i, core))), :i)
    @test 0 < v < vcore < min(L^3, sphere_vol(0.25))
    b = Box(center=c, size=(0.5, 0.5, 0.5))
    @test volume(geo(RVE(unitcell(), Phase(:i, Cut(b, Sphere(c, 0.2))))), :i) ≈ 0.125 - sphere_vol(0.2) rtol = 1e-6
    @test volume(geo(RVE(unitcell(), Phase(:i, Intersect(b, Box((0.5, 0.5, 0.5), (1, 1, 1)))))), :i) ≈ 0.25^3 rtol = 1e-9
    # dilation turns the sphere into a B-spline: volume accurate to ~1e-4
    @test volume(geo(RVE(unitcell(), Phase(:i, Scale(Sphere(c, 0.1), 2)))), :i) ≈ 8 * sphere_vol(0.1) rtol = 1e-3
    @test volume(geo(RVE(unitcell(), Phase(:i, Scale(Sphere(c, 0.1), (2, 1, 1))))), :i) ≈ 2 * sphere_vol(0.1) rtol = 1e-3
    # translating by a lattice vector changes nothing
    s = Sphere((0.9, 0.5, 0.5), 0.2)
    @test volume(geo(RVE(unitcell(), Phase(:i, Translate(s, (1, 0, -1))))), :i) ≈ sphere_vol(0.2) rtol = 1e-6
    @test volume(geo(RVE(unitcell(), Phase(:i, Translate(s, (0.3, 0, 0))))), :i) ≈ sphere_vol(0.2) rtol = 1e-6
    @test volume(geo(RVE(unitcell(), Phase(:i, Sphere(c, 0.2) ∪ Sphere((0.7, 0.5, 0.5), 0.2)))), :i) <
          2 * sphere_vol(0.2)
end

@testset "transform centres (positions)" begin
    sc = Scale(Sphere((0.6, 0.6, 0.6), 0.05), 2; center=(0.5, 0.5, 0.5))         # -> (0.7, 0.7, 0.7)
    ro = Rotate(Sphere((0.7, 0.5, 0.5), 0.05); angle=π / 2, axis=(0, 0, 1), center=(0.5, 0.5, 0.5))  # -> (0.5, 0.7, 0.5)
    generate(RVE(unitcell(), [Phase(:sc, sc), Phase(:ro, ro)]); geometry_only=true) do res
        com(n) = gmsh.model.occ.getCenterOfMass(3, only(res[n].entities))
        @test all(isapprox.(com(:sc), (0.7, 0.7, 0.7); atol=1e-4))   # B-spline: COM to ~1e-5
        @test all(isapprox.(com(:ro), (0.5, 0.7, 0.5); atol=1e-9))
    end
end

@testset "overlapping phases (G8)" begin
    c = (0.5, 0.5, 0.5)
    r1, r2 = 0.2, 0.3
    res = geo(RVE(unitcell(), [Phase(:coating, Sphere(c, r2)), Phase(:core, Sphere(c, r1))]))
    @test volume(res, :core) ≈ sphere_vol(r1) rtol = 1e-6
    @test volume(res, :coating) ≈ sphere_vol(r2) - sphere_vol(r1) rtol = 1e-6
    hidden = @test_logs (:warn, r"no volume") geo(RVE(unitcell(), [Phase(:core, Sphere(c, r1)),
                                                                     Phase(:coating, Sphere(c, r2))]))
    @test volume(hidden, :core) == 0 && hidden[:core].tag == 0
    @test volume(hidden, :coating) ≈ sphere_vol(r2) rtol = 1e-6
    @test_throws ErrorException geo(RVE(unitcell(), [Phase(:a, Sphere(c, r2)), Phase(:b, Sphere(c, r1))];
                                        overlaps=:error))
    # overlapping items of the same phase are fine in :error mode
    ok = geo(RVE(unitcell(), Phase(:a, [Sphere(c, 0.2), Sphere((0.6, 0.5, 0.5), 0.2)]); overlaps=:error))
    @test sphere_vol(0.2) < volume(ok, :a) < 2 * sphere_vol(0.2)
end

@testset "matrix name and no matrix (G9)" begin
    res = geo(RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), r)); matrix=:resin))
    @test phase_names(res) == ["resin", "p"]
    @test volume_fraction(res, "resin") ≈ 1 - sphere_vol(r) rtol = 1e-6
    pores = geo(RVE(unitcell(), Phase(:solid, Sphere((0.5, 0.5, 0.5), 0.6)); matrix=nothing))
    @test phase_names(pores) == ["solid"]
    @test volume_fraction(pores, :solid) < 1
    @test !("matrix" in group_names(pores))
    @test_throws KeyError volume(pores, :matrix)
end

@testset "non-unit shifted cell and mixed periodicity (G10)" begin
    cell = Cell((2.0, 1.0, 0.5); origin=(-3.0, 7.0, 0.25))
    res = geo(RVE(cell, Phase(:p, Sphere(cellpoint(cell, (1, 0, 1)), 0.2))))
    @test volume(res, :p) ≈ sphere_vol(0.2) rtol = 1e-6
    @test volume(res) ≈ 1.0
    mixed = Cell((1.0, 1.0, 1.0); periodic=(true, true, false))
    zs = geo(RVE(mixed, Phase(:p, Sphere((0.5, 0.5, 0.0), r))))
    @test volume(zs, :p) ≈ sphere_vol(r) / 2 rtol = 1e-6
    xs = geo(RVE(mixed, Phase(:p, Sphere((0.0, 0.5, 0.5), r))))
    @test volume(xs, :p) ≈ sphere_vol(r) rtol = 1e-6
    # large and small cells at moderate scales
    for k in (1e-2, 1e2)
        ck = Cell((k, k, k); origin=(k, 0, -k))
        rk = geo(RVE(ck, Phase(:p, Sphere(cellpoint(ck, (0, 0, 0)), 0.2k))))
        @test volume_fraction(rk, :p) ≈ sphere_vol(0.2) rtol = 1e-6
    end
end

@testset "errors and warnings (G11)" begin
    @test_logs (:warn, r"outside the cell") (:warn, r"no volume") geo(RVE(unitcell(periodic=false),
                                                                            Phase(:p, Sphere((3, 3, 3), 0.1))))
    @test_throws ErrorException geo(RVE(unitcell(), Phase(:p, Cut(Sphere((0.5, 0.5, 0.5), 0.1),
                                                                 Sphere((0.5, 0.5, 0.5), 0.2)))))
    @test_throws ErrorException geo(RVE(unitcell(), Phase(:p, Intersect(Sphere((0.2, 0.5, 0.5), 0.1),
                                                                       Sphere((0.8, 0.5, 0.5), 0.1)))))
end
