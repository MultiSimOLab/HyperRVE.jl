using HyperRVE, Test, Random

@testset "open_pores" begin
    dist = truncated_lognormal(4.0; cv=0.35, min=2.0, max=8.0)
    f = laguerre_foam(; n=16, diameter=dist, porosity=0.9, rng=1)
    a = open_pores(f; open_fraction=0.3, rng=3)
    r = a.report
    @test r.windows == round(Int, 0.3 * r.walls)
    @test r.open_fraction ≈ r.windows / r.walls
    @test 0 < r.open_area_fraction < 1
    @test r.eligible <= r.walls && length(a.windows) == r.windows == length(r.window_diameters)
    @test all(d -> d >= 2f.wall, r.window_diameters)
    @test r.percolates == (true, true, true)
    @test r.closed_pores <= 2 && r.open_pores + r.closed_pores == 16
    @test r.porosity ≈ 0.9 rtol = 1e-6
    @test 0.8 < r.open_porosity <= r.porosity
    @test length(a.fluid) + length(a.closed) == length(pores(f)) + r.windows
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
    # exact air volume (polyhedral pores, flat walls): open pores + the windows across the walls,
    # also with wide windows (narrowed so that both ends open into the pores)
    for win in (0.4, 0.8)
        aw = open_pores(f; open_fraction=0.3, window=win, rng=3)
        rw = generate(RVE(f.cell, Phase(:air, aw.fluid); matrix=:walls); geometry_only=true)
        V = prod(f.cell.size)
        expected = aw.report.open_porosity + sum(π * d^2 / 4 * f.wall for d in aw.report.window_diameters) / V
        @test volume_fraction(rw, :air) ≈ expected rtol = 1e-6
    end
    # :area does not depend on the units
    sel(ff, s) = open_pores(ff; open_fraction=0.2, select=:area, rng=s).report.window_diameters
    k = 1e-3
    fk = Foam{3}(Cell(k .* f.cell.size), laguerre(Cell(k .* f.cell.size); seeds=[k .* x .* f.cell.size for x in f.tessellation.seeds],
                 volumes=k^3 .* grain_volumes(f.tessellation)), f.porosity, k * f.wall, f.rounding, k .* f.diameters, f.scale)
    @test sel(fk, 5) ≈ k .* sel(f, 5) rtol = 1e-6
    @test sel(f, 5) != sel(f, 6)
    # one fluid network: the other percolating networks are closed
    a63 = open_pores(f; open_fraction=0.12, rng=63)
    @test a63.report.open_pores + a63.report.closed_pores == 16
    # a cell not periodic along z: the fluid is the largest network, no window through the free faces
    cz = Cell((14.0, 14.0, 14.0); periodic=(true, true, false))
    fz = laguerre_foam(cz; diameter=dist, porosity=0.9, rng=3)
    az = open_pores(fz; open_fraction=0.3, rng=1)
    @test az.report.open_pores > 0 && az.report.percolates[3] == false
    @test all(c -> 0 < c.base[3] < 14 && 0 < c.base[3] + c.axis[3] < 14, az.windows)
    cn = Cell((14.0, 14.0, 14.0); periodic=false)
    fn = laguerre_foam(cn; diameter=dist, porosity=0.9, rng=3)
    @test open_pores(fn; open_fraction=0.3, rng=1).report.open_pores > 0
    # errors
    @test_throws r"open_fraction" open_pores(f; open_fraction=1.5)
    @test_throws r"window must" open_pores(f; window=1.0)
    @test_throws r"select" open_pores(f; select=:big)
    @test_throws r"isolated" open_pores(f; isolated=:drop)
    @test_throws r"can take a window" open_pores(f; open_fraction=0.9, min_window=2.0)
    @test_throws r"only 3D" open_pores(laguerre_foam(Cell((20.0, 20.0)); diameter=3.0, porosity=0.85, rng=1))
end
