# Triply periodic minimal surfaces: kernel (pure Julia) — surfaces, densities and levels.

const TP = HyperRVE
@testset "TPMS surfaces" begin
    for (kind, type, level) in ((:gyroid, :sheet, 0.46), (:schwarz_p, :network, -0.7), (:diamond, :network, -0.3))
        S = TP._tpms_surface(kind, type, level, 1, 24)
        @test TP._tp_watertight(S, 0:6) && TP._tp_watertight(S, (0, 7:12...))
        # solid and complement fill the cell
        @test TP._tp_solid_volume(S) + TP._tp_signed_volume(S, (7:12...,)) - TP._tp_signed_volume(S, 0:0) ≈ 1 atol = 1e-12
        # the facets on opposite faces coincide after the translation
        for d in 1:3
            o = [i for i in 1:3 if i != d]
            tris(face) = sort([sort([round.(S.P[i][o]; digits=9) for i in t]) for (t, l) in zip(S.tri, S.lab) if l == face || l == face + 6])
            @test tris(2d - 1) == tris(2d)
        end
    end
end

@testset "TPMS densities and levels" begin
    for (kind, type, ρ, ref) in ((:gyroid, :network, 0.3, -0.6148), (:gyroid, :sheet, 0.3, 0.4629),
                                 (:schwarz_p, :network, 0.3, -0.6978), (:schwarz_p, :sheet, 0.3, 0.5240))
        c = tpms_level(kind; type, density=ρ)
        @test c ≈ ref rtol = 0.02
        @test tpms_density(kind; type, level=c) ≈ ρ atol = 2e-4
    end
    # symmetry of Schwarz P: level 0 halves the cell
    @test tpms_density(:schwarz_p; type=:network, level=0.0) ≈ 0.5 atol = 1e-3
    # two periods: same density
    @test tpms_density(:schwarz_p; type=:network, level=-0.7, periods=2, resolution=12) ≈
          tpms_density(:schwarz_p; type=:network, level=-0.7, resolution=24) atol = 5e-3
    # percolation of the network
    @test TP._tpms_percolates(:schwarz_p, :network, -0.7, 1, 24)
    @test !TP._tpms_percolates(:schwarz_p, :network, -1.1, 1, 24)
end

@testset "TPMS errors" begin
    @test_throws ArgumentError tpms_level(:schwarz_p; type=:network, density=0.15)        # below the pinch-off
    @test_throws ArgumentError tpms_density(:gyroid; type=:sheet, level=0.02)             # too thin for the grid
    @test_throws ArgumentError tpms_density(:lidinoid; level=0.3)
    @test_throws ArgumentError tpms_density(:gyroid; type=:shell, level=0.3)
    @test_throws ArgumentError tpms_density(:gyroid; type=:sheet, level=-0.3)
    @test_throws ArgumentError tpms_level(:gyroid; density=1.2)
    @test_throws ArgumentError tpms_level(:gyroid; density=0.3, periods=0)
    @test_throws ArgumentError tpms_level(:gyroid; density=0.3, resolution=4)
end
