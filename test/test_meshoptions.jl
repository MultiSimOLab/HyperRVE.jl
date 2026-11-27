# Element order, algorithms, optimisation and quality statistics.

elem_order() = begin
    types, _, _ = gmsh.model.mesh.getElements(3)
    only(unique(gmsh.model.mesh.getElementProperties(t)[3] for t in types))
end

"Mesh volume of the elements of a phase (order-aware)."
function phase_mesh_volume(res, name)
    v = 0.0
    for e in res[name].entities
        _, et, _ = gmsh.model.mesh.getElements(3, e)
        for t in et
            v += sum(gmsh.model.mesh.getElementQualities(t, "volume"))
        end
    end
    return v
end

@testset "order 2: periodic and valid" begin
    c = (0.5, 0.5, 0.5)
    cases = [RVE(unitcell(), Phase(:p, Sphere((0, 0, 0), 0.25))),
             RVE(unitcell(), Phase(:f, [Fibre(c, (1, 0, 0), 0.2), Fibre((0, 1, 1), (1, 0, 0), 0.2)])),
             RVE(Cell((2.0, 1.0, 0.5); origin=(1.0, -3.0, 0.25)), Phase(:p, Sphere((3.0, -2.5, 0.5), 0.2)))]
    for rve in cases
        n1 = generate(rve; h=minimum(rve.cell.size) / 5).periodicity
        generate(rve; h=minimum(rve.cell.size) / 5, mesh=MeshOptions(order=2, optimize_curved=:elastic)) do res
            @test res.stats.order == 2 && elem_order() == 2
            @test all(c -> c.ok && c.max_error < 1e-8, res.periodicity)
            @test all(i -> res.periodicity[i].npairs > n1[i].npairs, eachindex(n1))
            for (d, x) in enumerate(("x", "y", "z"))
                @test brute_periodic_match("$(x)min", "$(x)max", d, rve.cell.size[d])
            end
            @test res.stats.min_jacobian > 0 && res.stats.ninverted == 0
            @test mesh_volume() ≈ volume(res) rtol = 1e-2
        end
    end
end

@testset "curved vs straight order 2, inverted elements" begin
    rve = RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.3)))
    "Largest distance of the nodes of the sphere surface from the exact sphere."
    function deviation(r)
        dev = 0.0
        for v in r[:p].entities, (_, s) in gmsh.model.getBoundary([(3, v)], false, false, false)
            gmsh.model.getType(2, abs(s)) == "Sphere" || continue
            _, x, _ = gmsh.model.mesh.getNodes(2, abs(s), true)
            X = reshape(x, 3, :)
            dev = max(dev, maximum(abs.(sqrt.(vec(sum(abs2, X .- 0.5; dims=1))) .- 0.3)))
        end
        return dev
    end
    generate(rve; h=0.15, mesh=MeshOptions(order=2)) do r
        @test deviation(r) < 1e-9                       # mid-nodes on the sphere
    end
    generate(rve; h=0.15, mesh=MeshOptions(order=2, curved=false)) do r
        @test deviation(r) > 1e-3                       # straight-sided: mid-nodes inside
        @test r.stats.min_jacobian > 0.99
    end
    # coarse curved mesh: inverted elements are reported (and warned), :elastic fixes them
    coarse = RVE(unitcell(), Phase(:p, Sphere((0, 0, 0), 0.25)))
    bad = @test_logs (:warn, r"[1-9][0-9]* inverted") match_mode = :any generate(coarse; h=0.2,
                                                       mesh=MeshOptions(order=2, curvature=6, optimize=:default, optimize_curved=:none))
    @test bad.stats.ninverted > 0
    good = generate(coarse; h=0.2, mesh=MeshOptions(order=2))          # :auto -> netgen + elastic
    @test good.options.optimize_curved == :elastic
    @test good.stats.ninverted == 0 && all(c -> c.ok, good.periodicity)
    # the canonical cases with the default size have no inverted elements by default
    for g in (Sphere((0.5, 0.5, 0.5), 0.25), Sphere((0, 0, 0), 0.25), Sphere((0.5, 0.5, 0.5), 0.06))
        r = generate(RVE(unitcell(), Phase(:p, g)); mesh=MeshOptions(order=2))
        @test r.stats.ninverted == 0 && all(c -> c.ok, r.periodicity)
    end
    @test_throws ArgumentError generate(coarse; h=0.2, mesh=MeshOptions(order=2, optimize_curved=:optimize))
    np = generate(RVE(unitcell(periodic=false), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.25))); h=0.2,
                  mesh=MeshOptions(order=2, optimize_curved=:optimize))
    @test np.stats.order == 2
end

@testset "order 2 written file and scaled cell" begin
    k = 1e-3
    f = joinpath(OUTDIR, "order2.msh")
    res = generate(RVE(Cell((k, k, k); origin=(2k, 0, 0)), Phase(:p, Sphere((2k, 0, 0), 0.3k))); h=0.25k,
                   mesh=MeshOptions(order=2), output=f)
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f)
        @test length(gmsh.model.mesh.getNodes()[1]) == res.stats.nnodes
        x = reshape(gmsh.model.mesh.getNodes()[2], 3, :)
        @test all(isapprox.(extrema(x[1, :]), (2k, 3k); rtol=1e-9))
        linked = count(e -> gmsh.model.mesh.getPeriodicNodes(e...)[1] != e[2], gmsh.model.getEntities(2))
        @test linked > 0
        @test elem_order() == 2
    finally
        gmsh.finalize()
    end
    # the $Periodic block lists the mid-nodes as well (not only the vertices)
    lines = readlines(f)
    i = findfirst(==("\$Periodic"), lines)
    npairs = 0
    j = i + 2
    for _ in 1:parse(Int, lines[i+1])
        j += 2
        n = parse(Int, lines[j]); npairs += n; j += n + 1
    end
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f)
        nface = 0
        for (d, t) in gmsh.model.getEntities(2)
            m, nt, _, _ = gmsh.model.mesh.getPeriodicNodes(d, t, true)
            m != t && (nface += length(nt))
        end
        @test npairs >= nface > 0
    finally
        gmsh.finalize()
    end
end

@testset "algorithms and optimisation" begin
    rve = RVE(unitcell(), Phase(:p, Sphere((1.0, 0.5, 0.5), 0.25); refine=0.05))
    for alg in (:delaunay, :frontal, :hxt)
        r = generate(rve; h=0.15, mesh=MeshOptions(algorithm=alg))
        @test all(c -> c.ok, r.periodicity) && r.stats.min_quality > 0
    end
    for salg in (:default, :frontal, :delaunay, :meshadapt, :automatic)
        r = generate(rve; h=0.15, mesh=MeshOptions(surface_algorithm=salg))
        @test all(c -> c.ok, r.periodicity) && r.stats.min_quality > 0
    end
    qn = generate(rve; h=0.15, mesh=MeshOptions(optimize=:none)).stats.min_quality
    qg = generate(rve; h=0.15, mesh=MeshOptions(optimize=:netgen))
    @test qg.stats.min_quality >= 0.8 * qn && all(c -> c.ok, qg.periodicity)
    ho = generate(rve; h=0.15, mesh=MeshOptions(order=2, optimize_curved=:elastic))
    @test ho.stats.ninverted == 0 && all(c -> c.ok, ho.periodicity)
    # 2D cells use the surface algorithm
    r2 = generate(RVE(Cell((1.0, 1.0)), Phase(:p, Disk((0.5, 0.5), 0.2))); h=0.1, mesh=MeshOptions(surface_algorithm=:delaunay, order=2))
    @test r2.stats.order == 2 && all(c -> c.ok, r2.periodicity)
end

@testset "quality statistics and report" begin
    res = generate(RVE(unitcell(), [Phase(:a, Sphere((0.3, 0.5, 0.5), 0.15)), Phase(:b, Sphere((0.75, 0.5, 0.5), 0.15))]); h=0.15)
    s = res.stats
    @test sum(s.quality_histogram) == s.nelements
    @test minimum(last, s.phase_quality) ≈ s.min_quality
    @test Set(first.(s.phase_quality)) == Set(["matrix", "a", "b"])
    @test all(0 .<= s.worst_point .<= 1) && s.worst_phase in ("matrix", "a", "b")
    @test isnan(s.min_jacobian) && s.ninverted == 0
    io = IOBuffer()
    quality_report(io, res)
    out = String(take!(io))
    @test occursin("[0.9, 1.0]", out) && occursin("worst element", out) && occursin("matrix", out)
    io2 = IOBuffer(); quality_report(io2, generate(RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.2))); geometry_only=true))
    @test occursin("no mesh", String(take!(io2)))
end

@testset "options reach Gmsh and statistics are right (mutations)" begin
    rve = RVE(unitcell(), [Phase(:a, Sphere((0.3, 0.5, 0.5), 0.15); refine=0.05), Phase(:b, Sphere((0.75, 0.5, 0.5), 0.15))])
    generate(rve; h=0.15, mesh=MeshOptions(algorithm=:frontal, surface_algorithm=:meshadapt, optimize=:netgen)) do r
        @test gmsh.option.getNumber("Mesh.Algorithm3D") == 4
        @test gmsh.option.getNumber("Mesh.Algorithm") == 1
        @test gmsh.option.getNumber("Mesh.OptimizeNetgen") == 1
    end
    generate(rve; h=0.15, mesh=MeshOptions(optimize=:none)) do r
        @test gmsh.option.getNumber("Mesh.Optimize") == 0
    end
    for opts in (MeshOptions(), MeshOptions(order=2, curvature=6, optimize=:default, optimize_curved=:none))
        generate(RVE(unitcell(), [Phase(:a, Sphere((0, 0, 0), 0.25)), Phase(:b, Sphere((0.5, 0.5, 0.5), 0.2))]);
                 h=0.2, mesh=opts) do r
            s = r.stats
            tags = Int[]; own = String[]
            for p in r.phases, e in p.entities
                _, et, _ = gmsh.model.mesh.getElements(3, e)
                for t in et
                    append!(tags, t); append!(own, fill(p.name, length(t)))
                end
            end
            q = gmsh.model.mesh.getElementQualities(tags, "minSICN")
            hist = zeros(Int, 10)
            for x in q
                hist[clamp(floor(Int, 10x) + 1, 1, 10)] += 1
            end
            @test collect(s.quality_histogram) == hist
            for (name, pq) in s.phase_quality
                @test pq == minimum(q[own .== name])
            end
            if s.order == 2
                sj = gmsh.model.mesh.getElementQualities(tags, "minSJ")
                @test s.min_jacobian == minimum(sj)
                @test s.ninverted == count(<=(0), sj)
                w = s.ninverted > 0 ? argmin(sj) : argmin(q)
            else
                @test s.ninverted == count(<=(0), q)
                w = argmin(q)
            end
            @test s.worst_phase == own[w]
            _, nodes, _, _ = gmsh.model.mesh.getElement(tags[w])
            c = sum(gmsh.model.mesh.getNode(n)[1] for n in nodes) ./ length(nodes)
            @test all(isapprox.(s.worst_point, Tuple(c); atol=1e-12))
        end
    end
    # caller-owned session: the meshing options are restored by name
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        for (o, v) in (("Mesh.Algorithm3D", 4), ("Mesh.Algorithm", 5), ("Mesh.OptimizeNetgen", 1), ("Mesh.ElementOrder", 2),
                       ("Mesh.SecondOrderLinear", 1), ("Mesh.HighOrderOptimize", 2), ("Mesh.HighOrderSavePeriodic", 0))
            gmsh.option.setNumber(o, v)
        end
        generate(rve; h=0.2, mesh=MeshOptions(order=2))
        @test gmsh.option.getNumber("Mesh.Algorithm3D") == 4 && gmsh.option.getNumber("Mesh.Algorithm") == 5
        @test gmsh.option.getNumber("Mesh.OptimizeNetgen") == 1 && gmsh.option.getNumber("Mesh.ElementOrder") == 2
        @test gmsh.option.getNumber("Mesh.SecondOrderLinear") == 1 && gmsh.option.getNumber("Mesh.HighOrderOptimize") == 2
        @test gmsh.option.getNumber("Mesh.HighOrderSavePeriodic") == 0
    finally
        gmsh.finalize()
    end
end

@testset "interior B-spline bodies are meshed (default surface algorithm)" begin
    c = (0.5, 0.5, 0.5)
    for g in (Ellipsoid(c, (0.35, 0.12, 0.12)), Rotate(Ellipsoid(c, (0.35, 0.12, 0.12)); axis=(1, 1, 0), angle=π / 4),
              Scale(Sphere(c, 0.1), (3.5, 1.2, 1.2)))
        for per in (true, false)
            generate(RVE(unitcell(periodic=per), Phase(:p, g)); h=0.15) do r
                @test r.stats.nelements > 0
                @test mesh_volume() ≈ 1 rtol = 1e-9
            end
        end
    end
    sup = generate(supercell(RVE(unitcell(), Phase(:p, Ellipsoid(c, (0.3, 0.1, 0.1)))), (2, 2, 1)); h=0.2)
    @test sup.stats.nelements > 0 && all(x -> x.ok, sup.periodicity)
end

@testset "B-spline surfaces next to other phases are meshed, incomplete meshes are errors" begin
    c = (0.5, 0.5, 0.5)
    e = Ellipsoid(c, (0.3, 0.12, 0.12))
    for (rve, label) in ((RVE(unitcell(), [Phase(:a, Box((0, 0, 0.2), (1, 1, 0.6))), Phase(:b, e)]), "slab"),
                         (RVE(unitcell(periodic=false), [Phase(:a, Box((0, 0, 0.2), (1, 1, 0.6))), Phase(:b, e)]), "slab np"),
                         (RVE(unitcell(), [Phase(:a, Box((0, 0, 0), (1, 1, 1))), Phase(:b, e)]; overlaps=:paint), "full box"),
                         (RVE(unitcell(), Phase(:solid, Cut(Box((0, 0, 0), (1, 1, 1)), e)); matrix=nothing), "pore"))
        generate(rve; h=0.12) do r
            @test mesh_volume() ≈ volume(r) - (label == "pore" ? 4π * 0.3 * 0.12 * 0.12 / 3 : 0.0) rtol = 2e-2
            for p in r.phases
                p.volume > 0 && @test !isempty(p.entities)
            end
        end
    end
end

@testset "sphere cap cut next to its pole: overlapping facets are remeshed" begin
    # centre 0.01 from the face x = 0: Frontal-Delaunay overlaps the cap and its face disk on
    # this Mac (knife-edge: other platforms may mesh it at once; both must give a valid mesh)
    c = Cell((1.0, 1.0, 1.0))
    rve = RVE(c, Phase(:p, [Sphere((0.01, 0.5, 0.5), 0.2)]; refine=0.02))
    vm = Ref(0.0)
    res = generate(rve; h=0.1) do r
        vm[] = phase_mesh_volume(r, "p")
    end
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test vm[] ≈ 4π * 0.2^3 / 3 rtol = 1e-2
    @test res.stats.min_quality > 0.01 && res.stats.ninverted == 0
    # a scaled, shifted cell
    k = 1e-3
    cs = Cell(k .* (1.0, 1.0, 1.0); origin=(2.0, -1.0, 0.5))
    rs = generate(RVE(cs, Phase(:p, [Sphere(cs.origin .+ k .* (0.01, 0.5, 0.5), 0.2k)]; refine=0.02k)); h=0.1k)
    @test all(x -> x.ok, rs.periodicity) && rs.stats.ninverted == 0
    # the caller's logger and its old errors: kept, and not taken as errors of this mesh
    gmsh.initialize()
    try
        gmsh.logger.start()
        gmsh.option.setNumber("General.AbortOnError", 0)
        gmsh.logger.write("an old failure", "error")
        gmsh.option.setNumber("General.AbortOnError", 2)
        @test generate(rve; h=0.1).stats.ninverted == 0
        @test any(m -> occursin("an old failure", m), gmsh.logger.get())
        @test gmsh.option.getNumber("General.AbortOnError") == 2
    finally
        gmsh.logger.stop()
        gmsh.finalize()
    end
end
