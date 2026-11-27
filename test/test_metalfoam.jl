using HyperRVE, Test, Random

@testset "truncated_lognormal" begin
    d = truncated_lognormal(4.0; cv=0.35, min=2.0, max=8.0)
    rng = Random.Xoshiro(1)
    x = [d(rng) for _ in 1:200_000]
    @test sum(x) / length(x) ≈ 4.0 rtol = 5e-3          # the mean of the truncated law
    @test all(v -> 2.0 <= v <= 8.0, x)
    @test d(Random.Xoshiro(3)) == d(Random.Xoshiro(3))
    # strongly truncated (the mean far from the untruncated one)
    d2 = truncated_lognormal(2.5; cv=0.6, min=2.0, max=3.0)
    y = [d2(rng) for _ in 1:100_000]
    @test sum(y) / length(y) ≈ 2.5 rtol = 5e-3
    @test truncated_lognormal(4.0; cv=0, min=2.0, max=8.0)(rng) == 4.0
    @test occursin("truncated_lognormal(4.0", sprint(show, d))
    @test_throws ArgumentError truncated_lognormal(4.0; min=5.0, max=8.0)
    @test_throws ArgumentError truncated_lognormal(4.0; min=8.0, max=2.0)
    @test_throws ArgumentError truncated_lognormal(4.0; cv=-1, min=2.0, max=8.0)
    @test d(5) == d(Random.Xoshiro(5))                     # an integer seed
    @test truncated_lognormal(4.0; cv=1e-5, min=2.0, max=8.0)(rng) ≈ 4.0 rtol = 1e-3
    @test_throws r"finite" truncated_lognormal(4.0; min=2.0, max=Inf)
    @test_throws r"fraction" truncated_lognormal(1.00002; cv=0.5, min=1.0, max=1.00005)
end

@testset "laguerre_foam" begin
    dist = truncated_lognormal(4.0; cv=0.35, min=2.0, max=8.0)
    # n pores, cubic cell sized to hold them; polyhedral pores with the wall of the porosity
    f = laguerre_foam(; n=40, diameter=dist, porosity=0.9, rng=1)
    @test f isa Foam{3} && length(f.diameters) == 40
    @test prod(f.cell.size) ≈ sum(π .* f.diameters .^ 3 ./ 6) / 0.9 rtol = 1e-10
    st = pore_statistics(f)
    @test st.n == 40 && f.rounding == 0
    @test st.porosity ≈ 0.9 rtol = 1e-8
    @test abs(st.mean / st.target.mean - 1) < 0.05            # a uniform wall shrinks small pores more
    @test 2.0 * 0.9 < st.min && st.max < 8.0
    @test 12 < st.faces_per_pore < 16
    @test length(pores(f)) == 40
    # given cell and wall: the pores are rounded to reach the porosity
    fr = laguerre_foam(f.cell; diameter=dist, porosity=0.9, wall=0.1, rng=2)
    @test fr.wall == 0.1 && 0 < fr.rounding < 1
    sr = pore_statistics(fr)
    @test sr.porosity ≈ 0.9 rtol = 3e-3
    res = generate(RVE(fr.cell, Phase(:pores, pores(fr); void=true); matrix=:walls); h=1.2, retries=3)
    @test 1 - volume_fraction(res, :walls) ≈ sr.porosity rtol = 1e-4
    @test all(x -> x.ok, res.periodicity)
    # given rounding: the wall is solved
    fw = laguerre_foam(f.cell; diameter=dist, porosity=0.9, rounding=0.95, rng=2)
    @test fw.rounding == 0.95
    @test pore_statistics(fw).porosity ≈ 0.9 rtol = 3e-3
    # a vector of diameters, monodisperse, 2D, determinism
    fv = laguerre_foam(; n=3, diameter=[3.0, 4.0, 5.0], porosity=0.85, rng=1)
    @test sort(fv.diameters) ≈ [3.0, 4.0, 5.0]
    fm = laguerre_foam(Cell((10.0, 10.0, 10.0)); diameter=3.0, porosity=0.9, rng=1)
    @test all(==(fm.diameters[1]), fm.diameters) && abs(fm.diameters[1] / 3 - 1) < 0.1
    f2 = laguerre_foam(Cell((20.0, 20.0)); diameter=truncated_lognormal(3.0; min=1.5, max=6.0), porosity=0.85, rng=4)
    @test f2 isa Foam{2}
    @test pore_statistics(f2).porosity ≈ 0.85 rtol = 1e-8
    @test laguerre_foam(; n=10, diameter=dist, porosity=0.9, rng=7).tessellation.seeds ==
          laguerre_foam(; n=10, diameter=dist, porosity=0.9, rng=7).tessellation.seeds
    # the pore of each grain gets its target volume (grain volume × porosity)
    @test grain_volumes(f.tessellation) .* 0.9 ≈ π .* f.diameters .^ 3 ./ 6 rtol = 1e-6
    @test f.scale == 1.0 || abs(f.scale - 1) < 1e-12
    # the polyhedral foam meshes with the porosity
    resp = generate(RVE(f.cell, Phase(:pores, pores(f); void=true); matrix=:walls); h=1.5, retries=3)
    @test 1 - volume_fraction(resp, :walls) ≈ 0.9 rtol = 1e-6
    # a cell that does not match the diameters: scaled, with a warning
    fs = @test_logs (:warn, r"scaled by") match_mode = :any laguerre_foam(Cell((100.0, 100.0, 100.0)); diameter=[30.0, 40.0, 50.0], porosity=0.9, rng=1)
    @test fs.scale > 1.5 && fs.diameters ≈ sort([30.0, 40.0, 50.0]; rev=true) .* fs.scale
    # a wide law at low porosity: the smallest grains vanish, with a warning; statistics without them
    wide = truncated_lognormal(3.0; cv=0.8, min=0.5, max=10.0)
    fl = @test_logs (:warn, r"pores vanish") match_mode = :any laguerre_foam(; n=40, diameter=wide, porosity=0.4, rng=1)
    sl = pore_statistics(fl)
    @test sl.vanished > 0 && sl.n + sl.vanished == 40 && sl.min > 0
    # few pores: robust volumes (several seeds)
    @test all(s -> laguerre_foam(; n=4, diameter=dist, porosity=0.9, rng=s) isa Foam, 10:20)
    # another scale and origin; a cell not periodic along z
    k = 1e-3
    fk = laguerre_foam(Cell(k .* f.cell.size; origin=(5.0, -3.0, 100.0)); diameter=k .* f.diameters, porosity=0.9, rng=1)
    @test fk.wall ≈ k * laguerre_foam(f.cell; diameter=f.diameters, porosity=0.9, rng=1).wall rtol = 1e-6
    cz = Cell((14.0, 14.0, 14.0); periodic=(true, true, false))
    fz = laguerre_foam(cz; diameter=dist, porosity=0.9, rng=3)
    @test pore_statistics(fz).porosity ≈ 0.9 rtol = 1e-6
    # 2D with rounding; packing up to 0.8 in 2D
    f2r = laguerre_foam(Cell((20.0, 20.0)); diameter=truncated_lognormal(3.0; min=1.5, max=6.0), porosity=0.85, rounding=0.9, packing=0.7, rng=4)
    @test pore_statistics(f2r).porosity ≈ 0.85 rtol = 5e-3
    # errors
    @test_throws r"out of reach" laguerre_foam(f.cell; diameter=dist, porosity=0.9, wall=1.0, rng=1)
    @test_throws r"positive" laguerre_foam(f.cell; diameter=dist, porosity=0.9, wall=-0.1, rng=1)
    @test_throws r"positive and finite" laguerre_foam(; n=3, diameter=[3.0, -1.0, 2.0], porosity=0.9)
    @test_throws r"positive and finite" laguerre_foam(; n=3, diameter=_ -> NaN, porosity=0.9)
    @test_throws r"does not fit" laguerre_foam(Cell((10.0, 10.0, 10.0)); diameter=[12.0, 1.0], porosity=0.5)
    @test_throws r"wall or rounding" laguerre_foam(f.cell; diameter=dist, porosity=0.9, wall=0.1, rounding=0.5)
    @test_throws ArgumentError laguerre_foam(f.cell; diameter=dist, porosity=1.2)
    @test_throws r"fewer than 2 pores" laguerre_foam(Cell((3.0, 3.0, 3.0)); diameter=8.0, porosity=0.9)
    @test_throws ArgumentError laguerre_foam(; n=1, diameter=dist, porosity=0.9)
    @test occursin("Foam{3}(40 pores", sprint(show, f))
end
