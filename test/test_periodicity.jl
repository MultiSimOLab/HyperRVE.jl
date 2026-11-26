@testset "periodic node matching (N1, N3)" begin
    c = (0.5, 0.5, 0.5)
    cases = [
        "corner sphere" => RVE(unitcell(), Phase(:p, Sphere((0, 0, 0), 0.25); refine=InterfaceRefinement(level=2))),
        "fcc" => RVE(unitcell(), Phase(:p, [Sphere((0, 0.5, 0.5), 0.25), Sphere((0.5, 0, 0.5), 0.25),
                                            Sphere((0.5, 0.5, 0), 0.25), Sphere((1, 1, 1), 0.25)])),
        "fibres" => RVE(unitcell(), Phase(:f, [Fibre((0.5, 0.5, 0.5), (1, 0, 0), 0.25), Fibre((0, 1, 1), (1, 0, 0), 0.25)])),
        "boolean" => RVE(unitcell(), Phase(:i, Cut(Intersect(Sphere(c, 0.25), Box(center=c, size=(0.38, 0.38, 0.38))),
                                                    Fuse(Fibre(c, (1, 0, 0), 0.1), Fibre(c, (0, 1, 0), 0.1),
                                                         Fibre(c, (0, 0, 1), 0.1))))),
        "anisotropic shifted" => RVE(Cell((2.0, 1.0, 0.5); origin=(1.0, -3.0, 0.25)),
                                     Phase(:p, Sphere((3.0, -2.5, 0.5), 0.2))),
    ]
    for (name, rve) in cases
        @testset "$name" begin
            generate(rve; h=minimum(rve.cell.size) / 6) do res
                @test length(res.periodicity) == 3
                for chk in res.periodicity
                    @test chk.ok && chk.npairs > 0 && chk.max_error < 1e-8
                    d = chk.direction
                    x = "xyz"[d]
                    @test brute_periodic_match("$(x)min", "$(x)max", d, rve.cell.size[d])
                end
                @test res.stats.nelements > 0 && res.stats.min_quality > 0
                @test mesh_volume() ≈ volume(res) rtol = 1e-2
                @test 0 < res.stats.hmin <= res.stats.hmax
            end
        end
    end
end

@testset "Gmsh periodic node map (N2)" begin
    rve = RVE(unitcell(), Phase(:p, Sphere((1, 0.5, 0.5), 0.25)))
    generate(rve; h=0.2) do res
        g = only(x for x in res.groups if x.name == "xmax" && x.dim == 2)
        s = first(g.entities)
        master, nodes, _, _ = gmsh.model.mesh.getPeriodicNodes(2, s)
        xmin = only(x for x in res.groups if x.name == "xmin" && x.dim == 2)
        @test master in xmin.entities
        @test !isempty(nodes)
    end
end

@testset "mixed periodicity and no check (N4, N5)" begin
    rve = RVE(Cell((1.0, 1.0, 1.0); periodic=(true, false, true)), Phase(:p, Sphere((0.5, 0.0, 0.0), 0.25)))
    res = generate(rve; h=0.2)
    @test [c.direction for c in res.periodicity] == [1, 3]
    @test all(c -> c.ok, res.periodicity)
    @test isempty(generate(rve; h=0.2, check_periodicity=false).periodicity)
    np = generate(RVE(unitcell(periodic=false), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.2))); h=0.25)
    @test isempty(np.periodicity) && np.stats.nelements > 0
end

@testset "interface refinement" begin
    s = Sphere((0.5, 0.5, 0.5), 0.2)
    coarse = generate(RVE(unitcell(), Phase(:p, s)); h=0.2)
    fine = generate(RVE(unitcell(), Phase(:p, s; refine=InterfaceRefinement(h=0.03))); h=0.2)
    @test fine.stats.nelements > 2 * coarse.stats.nelements
    @test fine.stats.hmin < coarse.stats.hmin
end
