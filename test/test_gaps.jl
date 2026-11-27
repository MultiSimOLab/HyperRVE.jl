# Thin gap refinement (MeshOptions(gaps)): detection, size law, callback meshes.

const HG = HyperRVE
@testset "gap detection and law" begin
    balls = [((0.25, 0.25, 0.1), 0.1), ((0.25, 0.25, 0.302), 0.1)]
    gs = HG._find_gaps(balls, (0.0, 0.0, 0.0), (0.5, 0.5, 0.404), (true, true, true), 0.1, false)
    t = sort([g.t for g in gs])
    @test t[1] ≈ 0.002 && t[end] ≈ 0.002                       # the gap and its image across the z faces
    g = only(x for x in gs if abs(x.p[3] - 0.201) < 1e-9)
    @test g.R ≈ 0.05 && all(abs.(g.u) .≈ (0, 0, 1))
    @test length(gs) == 3                                         # the image gap at both z faces
    # contact, overlap and wide gaps do not count; faces of a free cell
    @test isempty(HG._find_gaps([((0.2, 0.2), 0.1), ((0.4, 0.2), 0.1)], (0.0, 0.0), (1.0, 1.0), (false, false), 0.1, false))
    @test isempty(HG._find_gaps([((0.2, 0.2), 0.1), ((0.5, 0.2), 0.1)], (0.0, 0.0), (1.0, 1.0), (false, false), 0.05, false))
    gf = HG._find_gaps([((0.103, 0.5), 0.1)], (0.0, 0.0), (1.0, 1.0), (false, false), 0.05, true)
    @test length(gf) == 1 && only(gf).t ≈ 0.003 && only(gf).face == 1
    # scaled cell: the same gaps scaled
    gs2 = HG._find_gaps([(b[1] .* 1e-6 .+ (2.0, 0.0, 0.0), b[2] * 1e-6) for b in balls], (2.0, 0.0, 0.0), (0.5e-6, 0.5e-6, 0.404e-6), (true, true, true), 0.1e-6, false)
    @test sort([x.t for x in gs2]) ≈ 1e-6 .* t
    # the law: t/n at the midpoint, growing away from it
    G = HG._GapI((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.002, 0.1, 2.0, 0.001)
    @test HG._gap_law(G, (0.0, 0.0, 0.0), 0.3) ≈ 0.001
    @test HG._gap_law(G, (0.02, 0.0, 0.0), 0.3) > HG._gap_law(G, (0.01, 0.0, 0.0), 0.3) > 0.001
    @test HG._gap_law(G, (0.0, 0.0, 0.05), 0.3) > HG._gap_law(G, (0.0, 0.0, 0.01), 0.3)
    # API
    @test MeshOptions(gaps=2).gaps == MeshOptions(gaps=true).gaps == MeshOptions(gaps=(elements=2,)).gaps == GapRefinement()
    @test MeshOptions().gaps === nothing && MeshOptions(gaps=false).gaps === nothing
    @test HG._with(MeshOptions(gaps=3); h=0.1).gaps == GapRefinement(elements=3)
    @test_throws ArgumentError MeshOptions(gaps=0)
    @test_throws ArgumentError MeshOptions(gaps=:yes)
    @test_throws ArgumentError GapRefinement(growth=0.9)
    @test occursin("elements = 2", sprint(show, GapRefinement()))
end

# elements whose barycentre is within t of a point, and the mean of their largest edges
function ligament(p, t)
    n, edge = 0, 0.0
    for (ty, tags, nodes) in zip(gmsh.model.mesh.getElements(3)...)
        np = gmsh.model.mesh.getElementProperties(ty)[4]
        for k in 1:np:length(nodes)
            X = [gmsh.model.mesh.getNode(nodes[k+j])[1] for j in 0:np-1]
            c = sum(X) ./ np
            sqrt(sum(abs2, c .- p)) < t || continue
            n += 1
            edge += maximum(sqrt(sum(abs2, X[a] .- X[b])) for a in 1:np for b in a+1:np)
        end
    end
    return n, n == 0 ? Inf : edge / n
end

@testset "gap meshes" begin
    r = 0.1; g = 0.04r
    cell = Cell((0.5, 0.5, 2 * (2r + g)))
    items = [Sphere((0.25, 0.25, r + g / 2), r), Sphere((0.25, 0.25, 3r + 1.5g), r)]
    p = (0.25, 0.25, 2r + g)
    n0 = Ref(0); n1 = Ref((0, 0.0))
    generate(RVE(cell, Phase(:p, items)); h=0.1, mesh=MeshOptions(optimize=:netgen)) do _
        n0[] = first(ligament(p, g))
    end
    res = generate(RVE(cell, Phase(:p, items)); h=0.1, mesh=MeshOptions(gaps=2, optimize=:netgen)) do _
        n1[] = ligament(p, g)
        @test brute_periodic_match("zmin", "zmax", 3, cell.size[3])
    end
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume(res, :p) ≈ 2 * 4π * r^3 / 3 rtol = 1e-6
    @test n0[] < 10 && n1[][1] >= 50 && n1[][2] <= g
    @test res.stats.min_quality > 0.1 && res.stats.nelements < 100_000
    res0 = res.stats.nelements
    k = 1e-6
    cs6 = Cell(k .* cell.size; origin=(3.0, -2.0, 5.0))
    its6 = [Sphere(cs6.origin .+ k .* (s.center), k * s.radius) for s in items]
    @test abs(generate(RVE(cs6, Phase(:p, its6)); h=0.1k, mesh=MeshOptions(gaps=2, optimize=:netgen)).stats.nelements - res0) < 0.1res0
    # 2D: far fewer elements than refining the whole phase
    c2 = Cell((1.0, 1.0))
    d = [Disk((0.5, 0.25), 0.1), Disk((0.5, 0.452), 0.1)]
    rg = generate(RVE(c2, Phase(:p, d)); h=0.05, mesh=MeshOptions(gaps=2))
    rr = generate(RVE(c2, Phase(:p, d; refine=0.001)); h=0.05)
    @test all(x -> x.ok, rg.periodicity) && rg.stats.nelements < rr.stats.nelements / 3
    # gap cut by a periodic face; free cell with faces
    cs = Cell((0.5, 0.5, 2 * (2r + g)))
    its = [Sphere((0.25, 0.25, 0.0), r), Sphere((0.25, 0.25, 2r + g), r)]
    res = generate(RVE(cs, Phase(:p, its)); h=0.1, mesh=MeshOptions(gaps=2))
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    cf = Cell((0.5, 0.5, 0.5); periodic=false)
    for faces in (true, false)
        res = generate(RVE(cf, Phase(:p, Sphere((0.25, 0.25, r + 0.003), r))); h=0.1, mesh=MeshOptions(gaps=GapRefinement(faces=faces)))
        @test res.stats.nelements > 0
    end
    # the callback is removed: remeshing the kept model gives the coarse mesh
    a = generate(RVE(cell, Phase(:p, items)); h=0.1).stats.nelements
    generate(RVE(cell, Phase(:p, items)); h=0.1, mesh=MeshOptions(gaps=2), keep_open=true)
    gmsh.model.mesh.clear()
    gmsh.model.mesh.generate(3)
    @test abs(length(gmsh.model.mesh.getElements(3)[2][1]) - a) < 0.2a
    close_gmsh()
    # far or non-finite points never throw; threads are restored
    @test HG._gap_size(HG._GAP_INDEX[], (1e20, -1e20, 0.0), 0.3) == 0.3
    @test HG._gap_size(HG._GAP_INDEX[], (NaN, 0.0, Inf), 0.3) == 0.3
    gmsh.initialize(); gmsh.option.setNumber("Mesh.MaxNumThreads3D", 4)
    generate(RVE(cell, Phase(:p, items)); h=0.1, mesh=MeshOptions(gaps=2), keep_open=true)
    @test gmsh.option.getNumber("Mesh.MaxNumThreads3D") == 4
    close_gmsh()
    # no matrix: no effect; frontal: a warning
    @test_logs (:warn, r"without a matrix") match_mode = :any generate(RVE(cell, Phase(:p, items); matrix=nothing); h=0.1, mesh=MeshOptions(gaps=2))
    # far items and a scaled, shifted cell give the same mesh
    far = [Translate(s, (0.0, 0.0, 3 * cell.size[3])) for s in items]
    @test abs(generate(RVE(cell, Phase(:p, far)); h=0.1, mesh=MeshOptions(gaps=2, optimize=:netgen)).stats.nelements - res0) < 0.1res0
    # no spheres: a warning
    @test_logs (:warn, r"no spheres") match_mode = :any generate(RVE(unitcell(), Phase(:b, Box((0.2, 0.2, 0.2), (0.3, 0.3, 0.3)))); h=0.2, mesh=MeshOptions(gaps=2))
end
