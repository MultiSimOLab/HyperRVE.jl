# Mesh control: interface refinement, periodic size field, phase sizes, bounds.

using RVEs: _image_shifts, classify_interfaces

q(res) = res.stats.min_quality
"Number of xmin nodes within `r` of (0, y, z)."
function near_xmin(y, z, r)
    _, x, _ = gmsh.model.mesh.getNodes()
    X = reshape(x, 3, :)
    return count(i -> abs(X[1, i] - X[1, argmin(X[1, :])]) < 1e-9 && hypot(X[2, i] - y, X[3, i] - z) < r, axes(X, 2))
end
"Mean longest edge of the elements of the volumes of a phase."
function mean_edge(res, name)
    e = Int[]
    for v in res[name].entities
        _, et, _ = gmsh.model.mesh.getElements(3, v)
        append!(e, reduce(vcat, et))
    end
    me = gmsh.model.mesh.getElementQualities(e, "maxEdge")
    return sum(me) / length(me)
end

@testset "image shifts (pure Julia)" begin
    lo, ext, per = (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (true, true, true)
    @test isempty(_image_shifts([(0.4, 0.4, 0.4, 0.6, 0.6, 0.6)], 0.1, lo, ext, per))
    @test _image_shifts([(0.6, 0.4, 0.4, 0.98, 0.6, 0.6)], 0.1, lo, ext, per) == [(1, 0, 0)]
    @test Set(_image_shifts([(0.6, 0.6, 0.4, 0.98, 0.98, 0.6)], 0.1, lo, ext, per)) == Set([(1, 0, 0), (0, 1, 0), (1, 1, 0)])
    @test isempty(_image_shifts([(0.6, 0.4, 0.4, 0.98, 0.6, 0.6)], 0.1, lo, ext, (false, true, true)))
    @test _image_shifts([(2.6, 0.4, 0.4, 2.98, 0.6, 0.6)], 0.1, (2.0, 0.0, 0.0), (1.0, 2.0, 2.0), per) == [(1, 0, 0)]
    @test _image_shifts([(0.02, 0.4, 0.4, 0.3, 0.6, 0.6)], 0.1, lo, ext, per) == [(-1, 0, 0)]
end

@testset "interfaces exclude cell faces and same-phase surfaces" begin
    generate(RVE(unitcell(), [Phase(:f, Fibre((0.5, 0.5, 0.5), (1, 0, 0), 0.2)),
                              Phase(:s, [Sphere((0.5, 0.15, 0.15), 0.1), Sphere((0.55, 0.15, 0.15), 0.1)])]);
             geometry_only=true) do res
        part = RVEs.Partition(3, phase_names(res), [p.tag for p in res.phases], [Int32.(p.entities) for p in res.phases],
                              Dict{Int32,Tuple{Int,Int}}())
        # rebuild piece ownership from the phases (items are not needed here)
        for (k, p) in enumerate(res.phases), e in p.entities
            part.piece_item[Int32(e)] = (k - 1, 1)
        end
        on = RVEs._boundary_entities(res.rve, RVEs.Frame(res.rve))
        ifs = classify_interfaces(part, on)
        tags = Set(Int(s.tag) for s in ifs)
        faces = Set(Int(t) for ((d, t), ps) in on.planes if d == 2)
        @test isempty(intersect(tags, faces))
        # fibre: lateral surface only (one surface); spheres: outer surface of the union only
        area = sum(gmsh.model.occ.getMass(2, t) for t in tags)
        @test area > 2π * 0.2 * 1.0 * 0.99
        @test all(s -> s.sides[1][1] != s.sides[2][1], ifs)
    end
end

@testset "periodic size field (E7)" begin
    e7(per; ps=true, cell=unitcell(periodic=per)) =
        generate(RVE(cell, Phase(:p, Sphere(cellpoint(cell, (0.78, 0.5, 0.5)), 0.2 * cell.size[1]); refine=0.02 * cell.size[1]));
                 h=0.15 * cell.size[1], mesh=MeshOptions(periodic_sizes=ps))
    ref = e7(false)
    per = e7(true)
    @test q(per) > 0.05
    @test q(per) >= 0.3 * q(ref)
    @test all(c -> c.ok, per.periodicity)
    # the refinement reaches xmin through the periodic images
    generate(RVE(unitcell(), Phase(:p, Sphere((0.78, 0.5, 0.5), 0.2); refine=0.02)); h=0.15) do r
        @test near_xmin(0.5, 0.5, 0.15) >= 20
    end
    generate(RVE(unitcell(), Phase(:p, Sphere((0.78, 0.5, 0.5), 0.2))); h=0.15) do r
        @test near_xmin(0.5, 0.5, 0.15) < 10
    end
    # without the periodic field xmin is not refined (keeps this test meaningful)
    generate(RVE(unitcell(), Phase(:p, Sphere((0.78, 0.5, 0.5), 0.2); refine=0.02)); h=0.15,
             mesh=MeshOptions(periodic_sizes=false)) do r
        @test near_xmin(0.5, 0.5, 0.15) < 10
        @test q(r) < q(per)
    end
    # lattice-shift invariance: same medium with the inclusion in the middle of the cell
    mid = generate(RVE(unitcell(), Phase(:p, Sphere((0.28, 0.5, 0.5), 0.2); refine=0.02)); h=0.15)
    @test abs(mid.stats.nelements - per.stats.nelements) <= 0.2 * per.stats.nelements
    # corner: refined sphere close to the xmax/ymax/zmax corner without touching it
    generate(RVE(unitcell(), Phase(:p, Sphere((0.83, 0.83, 0.83), 0.12); refine=0.02)); h=0.15) do r
        @test q(r) > 0.05
        @test all(c -> c.ok && c.npairs > 0, r.periodicity)
        for (d, x) in enumerate(("x", "y", "z"))
            @test brute_periodic_match("$(x)min", "$(x)max", d, 1.0)
        end
    end
    # shifted, non-cubic, tiled domain, and other scales
    c = Cell((1.0, 0.8, 1.2); origin=(3.0, -1.0, 0.5))
    tiled = generate(RVE(c, Phase(:p, Sphere(cellpoint(c, (0.8, 0.5, 0.5)), 0.15); refine=0.02); tiling=(2, 1, 1)); h=0.15)
    @test q(tiled) > 0.05
    @test all(x -> x.ok, tiled.periodicity)
    big = e7(true; cell=Cell((1e3, 1e3, 1e3)))
    @test abs(big.stats.nelements - per.stats.nelements) <= 0.01 * per.stats.nelements
    # semi-periodic: nothing is mirrored across the non-periodic direction
    semi = generate(RVE(unitcell(periodic=(false, true, true)), Phase(:p, Sphere((0.78, 0.5, 0.5), 0.2); refine=0.02)); h=0.15)
    @test q(semi) > 0.05
end

@testset "sizes: refine, growth, phase h, hmin, hmax" begin
    s = Sphere((0.5, 0.5, 0.5), 0.25)
    plain = generate(RVE(unitcell(), Phase(:p, s)); h=0.15)
    fine = generate(RVE(unitcell(), Phase(:p, s; refine=0.03)); h=0.15)
    @test fine.stats.hmin < 0.05 && fine.stats.nelements > plain.stats.nelements
    g12 = generate(RVE(unitcell(), Phase(:p, s; refine=(h=0.03, growth=1.2))); h=0.15).stats.nelements
    g2 = generate(RVE(unitcell(), Phase(:p, s; refine=(h=0.03, growth=2.0))); h=0.15).stats.nelements
    @test g12 > g2
    generate(RVE(unitcell(), Phase(:p, s; h=0.04)); h=0.15) do r
        @test mean_edge(r, :p) < 0.07
        @test mean_edge(r, :matrix) > mean_edge(r, :p)
    end
    generate(RVE(unitcell(), Phase(:p, s; h=0.2)); h=0.05) do r       # coarse inclusion, fine matrix
        @test mean_edge(r, :p) > mean_edge(r, :matrix)
    end
    bounded = generate(RVE(unitcell(), Phase(:p, s; refine=0.01)); h=0.15, mesh=MeshOptions(hmin=0.03))
    @test bounded.stats.hmin >= 0.5 * 0.03
    capped = generate(RVE(unitcell(), Phase(:p, s)); mesh=MeshOptions(h=0.3, hmax=0.1))
    @test capped.stats.hmax <= 2.5 * 0.1      # longest tet edge vs target size
    @test_logs (:warn, r"no effect") match_mode = :any generate(RVE(unitcell(), Phase(:p, s; refine=0.5)); h=0.15)
    # level rule of the 2024 API still works
    lv = generate(RVE(unitcell(), Phase(:p, s; refine=InterfaceRefinement(level=2))); h=0.15)
    @test lv.stats.nelements > plain.stats.nelements
end
