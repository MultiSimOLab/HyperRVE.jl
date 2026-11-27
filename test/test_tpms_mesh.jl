# TPMS phases (experimental): discrete surface meshed by Gmsh, groups, periodicity, $Periodic.

# node pairs of the $Periodic section of a .msh file: every pair is a translate by one cell edge
function periodic_pairs_ok(path, L)
    lines = readlines(path)
    i = findfirst(==("\$Nodes"), lines); j = findfirst(==("\$EndNodes"), lines)
    X = Dict{Int,NTuple{3,Float64}}()
    k = i + 2
    while k < j
        n = parse(Int, split(lines[k])[4])
        for m in 1:n
            X[parse(Int, lines[k+m])] = Tuple(parse.(Float64, split(lines[k+n+m])[1:3]))
        end
        k += 2n + 1
    end
    i = findfirst(==("\$Periodic"), lines)
    i === nothing && return (false, 0)
    nl = parse(Int, lines[i+1]); k = i + 2; ok = true; tot = 0
    slaves = Set{Int}()
    for _ in 1:nl
        k += 2
        np = parse(Int, lines[k]); k += 1
        for m in 0:np-1
            a, b = parse.(Int, split(lines[k+m]))
            δ = X[a] .- X[b]
            ok &= count(x -> abs(x - L) < 1e-9 * L, δ) == 1 && count(x -> abs(x) < 1e-9 * L, δ) == 2
            ok &= !(a in slaves)                 # each node is a slave once
            push!(slaves, a); tot += 1
        end
        k += np
    end
    return ok, tot
end

@testset "TPMS meshes" begin
    c = unitcell()
    g = tpms(:gyroid; type=:sheet, density=0.3)
    @test g isa TPMS && g.level ≈ tpms_level(:gyroid; type=:sheet, density=0.3)
    out = joinpath(OUTDIR, "gyroid.msh")
    res = generate(RVE(c, Phase(:gyroid, g)); output=out)
    @test volume_fraction(res, :gyroid) ≈ 0.3 atol = 2e-3
    @test volume_fraction(res, :gyroid) + volume_fraction(res, :matrix) ≈ 1 atol = 1e-12
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test res.stats.min_quality > 0.005
    @test "interface_matrix_gyroid" in interface_names(res)
    @test all(n -> n in group_names(res), ("xmin", "zmax", "corners", "boundary", "edge_xmin_ymin", "corner_xmax_ymax_zmax"))
    ok, n = periodic_pairs_ok(out, 1.0)
    @test ok && n > 0
    # no matrix, two periods, a small cell away from the origin
    k = 1e-3
    cs = Cell((k, k, k); origin=(2k, -k, 0.5k))
    p = tpms(:schwarz_p; type=:network, density=0.35, periods=2, resolution=16)
    out = joinpath(OUTDIR, "p.msh")
    res = generate(RVE(cs, Phase(:p, p); matrix=nothing); output=out)
    @test volume_fraction(res, :p) ≈ 0.35 atol = 3e-3
    @test volume(res, :p) ≈ 0.35 * k^3 rtol = 1e-2
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test "interface_p_void" in interface_names(res)
    ok, n = periodic_pairs_ok(out, k)
    @test ok && n > 0
    # geometry only: kernel fractions
    res = generate(RVE(c, Phase(:gyroid, g)); geometry_only=true)
    @test volume_fraction(res, :gyroid) ≈ 0.3 atol = 2e-3
end

@testset "TPMS errors" begin
    c = unitcell()
    g = tpms(:schwarz_p; type=:network, level=-0.7)
    @test_throws ArgumentError generate(RVE(c, [Phase(:p, g), Phase(:s, Sphere((0.5, 0.5, 0.5), 0.1))]))
    @test_throws ArgumentError generate(RVE(c, Phase(:p, [g, Sphere((0.5, 0.5, 0.5), 0.1)])))
    @test_throws ArgumentError generate(RVE(Cell((1.0, 1.0, 2.0)), Phase(:p, g)))
    @test_throws ArgumentError generate(RVE(Cell((1.0, 1.0, 1.0); periodic=(true, true, false)), Phase(:p, g)))
    @test_throws ArgumentError generate(RVE(c, Phase(:p, g); tiling=2))
    @test_throws ArgumentError generate(RVE(c, Phase(:p, g; refine=0.02)))
    @test_throws ArgumentError generate(RVE(c, Phase(:p, g)); mesh=MeshOptions(order=2))
    @test_throws ArgumentError tpms(:gyroid; type=:sheet)
    @test_throws ArgumentError tpms(:gyroid; type=:sheet, density=0.3, level=0.4)
end
