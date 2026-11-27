using HyperRVE, Test, Random

@testset "open_pores" begin
    dist = truncated_lognormal(4.0; cv=0.35, min=2.0, max=8.0)
    f = laguerre_foam(; n=16, diameter=dist, porosity=0.9, rng=1)
    a = open_pores(f; open_fraction=0.3, rng=3)
    r = a.report
    @test r.windows == round(Int, 0.3 * r.walls) || r.closed_pores > 0
    @test r.open_fraction ≈ round(Int, 0.3 * r.walls) / r.walls
    @test 0 < r.open_area_fraction < 1
    @test r.eligible <= r.walls && length(a.windows) == r.windows == length(r.window_diameters)
    @test all(d -> d >= 2f.wall, r.window_diameters)
    @test r.percolates == (true, true, true) && r.components == 1 && r.closed_pores == 0
    @test r.open_porosity ≈ 0.9 rtol = 1e-6
    @test r.porosity ≈ 0.9 rtol = 1e-6
    @test isempty(a.closed) && length(a.fluid) == length(pores(f)) + r.windows
    # determinism; more windows, fewer components
    @test open_pores(f; open_fraction=0.3, rng=3).report.window_diameters == r.window_diameters
    a1 = open_pores(f; open_fraction=0.1, rng=2)
    @test a1.report.components > 1 && a1.report.closed_pores > 0 && !isempty(a1.closed)
    @test a1.report.open_porosity < 0.9
    @test length(a1.fluid) + length(a1.closed) == length(pores(f)) + a1.report.windows
    # isolated pores kept on request
    ak = open_pores(f; open_fraction=0.1, rng=2, isolated=:keep)
    @test isempty(ak.closed)
    @test ak.report.open_porosity ≈ 0.9 rtol = 1e-6
    # no percolation: a warning
    @test_logs (:warn, r"percolate") match_mode = :any open_pores(f; open_fraction=0.02, rng=1)
    # select by area: larger open area for the same number
    @test open_pores(f; open_fraction=0.2, rng=4, select=:area).report.open_area_fraction >
          open_pores(f; open_fraction=0.2, rng=4).report.open_area_fraction
    # the fluid meshes: with the walls and alone; its volume is the open porosity plus the windows
    res = generate(RVE(f.cell, [Phase(:air, a1.fluid), Phase(:closed, a1.closed; void=true)]; matrix=:walls); h=1.5, retries=3)
    @test all(x -> x.ok, res.periodicity)
    @test volume_fraction(res, :air) > a1.report.open_porosity
    @test volume_fraction(res, :air) + volume_fraction(res, :walls) ≈ 1 - (a1.report.porosity - a1.report.open_porosity) rtol = 1e-3
    fr = laguerre_foam(f.cell; diameter=dist, porosity=0.9, wall=0.1, rng=1)
    ar = open_pores(fr; open_fraction=0.3, rng=5)
    res = generate(RVE(fr.cell, Phase(:air, ar.fluid); matrix=nothing); h=1.5, retries=3)
    @test "interface_air_void" in interface_names(res)
    @test volume_fraction(res, :air) > ar.report.open_porosity
    # errors
    @test_throws r"open_fraction" open_pores(f; open_fraction=1.5)
    @test_throws r"window must" open_pores(f; window=1.0)
    @test_throws r"select" open_pores(f; select=:big)
    @test_throws r"isolated" open_pores(f; isolated=:drop)
    @test_throws r"can take a window" open_pores(f; open_fraction=0.9, min_window=2.0)
    @test_throws r"only 3D" open_pores(laguerre_foam(Cell((20.0, 20.0)); diameter=3.0, porosity=0.85, rng=1))
end
