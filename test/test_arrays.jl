# Regular particle arrays and laminates.

@testset "particle arrays" begin
    c = unitcell()
    for (k, n) in ((:sc, 1), (:bcc, 2), (:fcc, 4))
        p = particle_array(c, k; fraction=0.3)
        @test length(p) == n && volume_fraction(p, c) ≈ 0.3
        @test volume_fraction(geo(RVE(c, Phase(:p, p))), :p) ≈ 0.3 rtol = 1e-8
    end
    @test particle_array(c, :sc; radius=0.2)[1].center == (0.5, 0.5, 0.5)
    # touching limits: sc 0.5236, bcc 0.6802, fcc 0.7405
    @test_throws ArgumentError particle_array(c, :sc; fraction=0.53)
    @test_throws ArgumentError particle_array(c, :bcc; fraction=0.69)
    @test length(particle_array(c, :fcc; fraction=0.74)) == 4
    @test_throws ArgumentError particle_array(c, :fcc; fraction=0.75)
    # any cell size and origin
    cb = Cell((2.0, 1.0, 1.0); origin=(1.0, 1.0, 1.0))
    @test all(s -> s.center[1] in (1.0, 2.0), particle_array(cb, :bcc; radius=0.2))
    # 2D
    c2 = Cell((1.0, sqrt(3)))
    res = generate(RVE(c2, Phase(:f, particle_array(c2, :hex; fraction=0.6))); h=0.05)
    @test volume_fraction(res, :f) ≈ 0.6 rtol = 1e-8
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test length(particle_array(Cell((1.0, 1.0)), :square; fraction=0.5)) == 1
    @test_throws ArgumentError particle_array(Cell((1.0, 1.0)), :hex; fraction=0.5)          # proportions
    @test_throws ArgumentError particle_array(c, :hex; fraction=0.3)                          # 2D kind in 3D
    @test_throws ArgumentError particle_array(c, :sc)
    @test_throws ArgumentError particle_array(c, :sc; radius=0.1, fraction=0.1)
end

@testset "laminates" begin
    c = unitcell()
    ph = laminate(c; axis=3, fractions=(0.3, 0.7))                                         # first layer centred on the face
    @test [p.name for p in ph] == ["layer_1", "layer_2"]
    res = generate(RVE(c, ph; matrix=nothing); h=0.1)
    @test [p.fraction for p in res.phases] ≈ [0.3, 0.7] rtol = 1e-8
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test interface_area(res, :layer_1, :layer_2) ≈ 2.0 rtol = 1e-8                     # two planes per period
    ph = laminate(Cell((2.0, 1.0)); axis=1, fractions=[0.2, 0.3, 0.5], names=(:a, :b, :c))
    g = geo(RVE(Cell((2.0, 1.0)), ph; matrix=nothing))
    @test [volume(g, n) for n in (:a, :b, :c)] ≈ [0.4, 0.6, 1.0] rtol = 1e-10
    @test_throws ArgumentError laminate(c; fractions=(0.3, 0.6))
    @test_throws ArgumentError laminate(c; fractions=(1.0,))
    @test_throws ArgumentError laminate(c; axis=4, fractions=(0.5, 0.5))
    @test_throws ArgumentError laminate(c; fractions=(0.5, 0.5), names=(:a,))
end

@testset "arrays and laminates: review fixes" begin
    # a laminate along a free axis fills the cell (layers crossing the face are split)
    cf = Cell((1.0, 1.0, 1.0); periodic=(true, true, false))
    g = geo(RVE(cf, laminate(cf; axis=3, fractions=(0.3, 0.7), offset=0.5); matrix=nothing))
    @test [p.fraction for p in g.phases] ≈ [0.3, 0.7] rtol = 1e-10
    @test !("interface_layer_1_void" in interface_names(g))
    # repeated names form one phase
    ph = laminate(unitcell(); fractions=(0.25, 0.25, 0.25, 0.25), names=(:a, :b, :a, :b))
    @test [p.name for p in ph] == ["a", "b"] && length(ph[1].items) >= 2
    g = geo(RVE(unitcell(), ph; matrix=nothing))
    @test volume(g, :a) ≈ 0.5 rtol = 1e-10
    # default offset: every interface plane is inside the cell
    res = generate(RVE(unitcell(), laminate(unitcell(); fractions=fill(0.1, 10)); matrix=nothing); h=0.2)
    @test length(interface_names(res)) == 10
    @test_throws ArgumentError laminate(unitcell(); fractions=(0.5, 0.5), names="ab")
    @test_throws ArgumentError laminate(unitcell(); fractions=("a", "b"))
    @test_throws ArgumentError laminate(unitcell())
    # particle arrays need a periodic cell when particles sit on faces
    @test_throws ArgumentError particle_array(Cell((1.0, 1.0, 1.0); periodic=false), :bcc; fraction=0.2)
    @test length(particle_array(Cell((1.0, 1.0, 1.0); periodic=false), :sc; radius=0.45)) == 1
    @test_throws ArgumentError particle_array(Cell((1.0, 1.0, 1.0); periodic=false), :sc; radius=0.6)
end
