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
    bad = @test_logs (:warn, r"inverted") match_mode = :any generate(coarse; h=0.2, mesh=MeshOptions(order=2))
    @test bad.stats.ninverted > 0
    good = generate(coarse; h=0.2, mesh=MeshOptions(order=2, optimize_curved=:elastic))
    @test good.stats.ninverted == 0 && all(c -> c.ok, good.periodicity)
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
end

@testset "algorithms and optimisation" begin
    rve = RVE(unitcell(), Phase(:p, Sphere((1.0, 0.5, 0.5), 0.25); refine=0.05))
    for alg in (:delaunay, :frontal, :hxt)
        r = generate(rve; h=0.15, mesh=MeshOptions(algorithm=alg))
        @test all(c -> c.ok, r.periodicity) && r.stats.min_quality > 0
    end
    for salg in (:frontal, :delaunay, :meshadapt, :automatic)
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
    @test occursin("[0.9, 1.0)", out) && occursin("worst element", out) && occursin("matrix", out)
    io2 = IOBuffer(); quality_report(io2, generate(RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.2))); geometry_only=true))
    @test occursin("no mesh", String(take!(io2)))
end
