# Capsule primitive and fused phases (kept cheap: small bodies, coarse meshes).

capsule_vol(r, L) = π * r^2 * L + 4π * r^3 / 3

"Manual BCC strut lattice (cylinders + node spheres), struts listed once."
function bcc_items(r)
    c = (0.5, 0.5, 0.5)
    items = Geometry{3}[Sphere(c, r), Sphere((0.0, 0.0, 0.0), r)]
    for k in Iterators.product((0, 1), (0, 1), (0, 1))
        push!(items, Cylinder(c, Float64.(k) .- c, r))
    end
    return items
end

@testset "Capsule (K1)" begin
    @test Capsule(center=(0, 0, 0), direction=(0, 0, 2), length=1, radius=0.1) == Capsule((0, 0, -0.5), (0, 0, 0.5), 0.1)
    @test_throws ArgumentError Capsule((0, 0, 0), (0, 0, 0), 0.1)
    @test_throws ArgumentError Capsule((0, 0, 0), (0, 0, 1), -0.1)
    @test RVEs.refpoint(Capsule((0, 0, 0), (0, 0, 1), 0.1)) == (0.0, 0.0, 0.5)
    for (p, q) in (((0.3, 0.5, 0.5), (0.7, 0.5, 0.5)), ((0.3, 0.3, 0.3), (0.7, 0.6, 0.5)), ((0.5, 0.5, 0.7), (0.5, 0.5, 0.3)))
        L = sqrt(sum(abs2, q .- p))
        generate(RVE(unitcell(), Phase(:c, Capsule(p, q, 0.1))); geometry_only=true) do res
            @test volume(res, :c) ≈ capsule_vol(0.1, L) rtol = 1e-9
            types = Set(gmsh.model.getType(2, t) for (_, t) in gmsh.model.getEntities(2))
            @test !("BSpline surface" in types)            # analytic surfaces only
        end
    end
    res = generate(RVE(unitcell(), Phase(:c, Capsule((0.8, 0.9, 0.5), (1.2, 1.1, 0.5), 0.08))); h=0.12)
    @test volume(res, :c) ≈ capsule_vol(0.08, sqrt(0.16 + 0.04)) rtol = 1e-6
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    k = 1e-3
    ck = Cell((k, k, k); origin=(2k, 0, 0))
    rk = geo(RVE(ck, Phase(:c, Capsule(cellpoint(ck, (0.3, 0.5, 0.5)), cellpoint(ck, (0.7, 0.5, 0.5)), 0.1k))))
    @test volume_fraction(rk, :c) ≈ capsule_vol(0.1, 0.4) rtol = 1e-9
    @test_logs (:warn, r"Capsule longer") RVE(unitcell(), Phase(:c, Capsule((-0.2, 0.5, 0.5), (1.2, 0.5, 0.5), 0.1)))
end

@testset "fused strut lattice (K2-K4)" begin
    r = 0.08
    un = geo(RVE(unitcell(), Phase(:lat, bcc_items(r)); matrix=nothing))
    @test length(un[:lat].entities) > 10
    local vf
    generate(RVE(unitcell(), Phase(:lat, bcc_items(r); fuse=true); matrix=nothing); h=0.1) do res
        vf = volume(res, :lat)
        @test vf ≈ volume(un, :lat) rtol = 1e-6
        @test length(res[:lat].entities) == 1
        for (_, s) in gmsh.model.getEntities(2)                     # no internal surfaces
            @test length(gmsh.model.getAdjacencies(2, s)[1]) == 1
        end
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
        for (d, x) in enumerate(("x", "y", "z"))
            @test brute_periodic_match("$(x)min", "$(x)max", d, 1.0)
        end
    end
    @test occursin("fused", sprint(show, Phase(:lat, bcc_items(r); fuse=true)))
    gi = generate(RVE(unitcell(), Phase(:lat, bcc_items(r); fuse=true)); h=0.1)
    @test interface_names(gi) == ["interface_matrix_lat"] && all(x -> x.ok, gi.periodicity)
end

@testset "fused supercells: unit path and flat path (K5, K6)" begin
    spheres = [Sphere((0.0, 0.0, 0.0), 0.3), Sphere((0.5, 0.5, 0.5), 0.2), Sphere((0.6, 0.5, 0.5), 0.2)]
    one = geo(RVE(unitcell(), Phase(:p, spheres; fuse=true); matrix=nothing))
    for n in ((2, 2, 2), (3, 1, 2))
        generate(supercell(RVE(unitcell(), Phase(:p, spheres; fuse=true); matrix=nothing), n); geometry_only=true) do s
            @test volume(s, :p) ≈ prod(n) * volume(one, :p) rtol = 1e-6
            for (_, f) in gmsh.model.getEntities(2)                 # tiles fused: no internal planes
                @test length(gmsh.model.getAdjacencies(2, f)[1]) == 1
            end
        end
    end
    # a tiled non-periodic direction must use the flat path (items sticking out of a cell
    # appear whole in the next copy); the unit path would clip them
    cell = Cell((1.0, 1.0, 1.0); periodic=(false, true, true))
    f = geo(RVE(cell, Phase(:p, [Sphere((1.0, 0.5, 0.5), 0.2), Sphere((0.8, 0.5, 0.5), 0.15)]; fuse=true); tiling=(2, 1, 1)))
    @test volume(f, :p) ≈ 0.0683787 rtol = 1e-4
    # 2D supercell
    d2 = [Disk((0.4, 0.5), 0.15), Disk((0.6, 0.5), 0.15), Disk((0.0, 0.0), 0.2)]
    a1 = geo(RVE(Cell((1.0, 1.0)), Phase(:p, d2; fuse=true)))
    a2 = geo(supercell(RVE(Cell((1.0, 1.0)), Phase(:p, d2; fuse=true)), (2, 3)))
    @test volume(a2, :p) ≈ 6 * volume(a1, :p) rtol = 1e-9
end

@testset "duplicates, fibres, refinement, errors (K7, K8)" begin
    # exact duplicate placements are dropped before the fuse (also with -0.0 coordinates)
    dup = [Sphere((0.5, 0.5, 0.5), 0.2), Sphere((0.5, 0.5, 0.5), 0.2), Cylinder((0.2, 0.5, 0.5), (0.6, 0.0, 0.0), 0.05),
           Cylinder((0.8, 0.5, 0.5), (-0.6, -0.0, 0.0), 0.05)]
    @test RVEs._placement_key(dup[3], (0.0, 0.0, 0.0)) == RVEs._placement_key(dup[4], (0.0, 0.0, 0.0))
    @test RVEs._placement_key(Sphere((-0.0, 0.5, 0.5), 0.1), (0.0, 0.0, 0.0)) == RVEs._placement_key(Sphere((0.0, 0.5, 0.5), 0.1), (0.0, 0.0, 0.0))
    single = geo(RVE(unitcell(), Phase(:p, dup[[1, 3]]; fuse=true)))
    both = geo(RVE(unitcell(), Phase(:p, dup; fuse=true)))
    @test volume(both, :p) ≈ volume(single, :p) rtol = 1e-9
    # fused aligned fibres crossing each other
    fib = geo(RVE(unitcell(), Phase(:f, [Fibre((0.5, 0.5, 0.5), (1, 0, 0), 0.1), Fibre((0.5, 0.5, 0.5), (0, 1, 0), 0.1)]; fuse=true)))
    @test π * 0.01 < volume(fib, :f) < 2π * 0.01 && length(fib[:f].entities) == 1
    # level refinement of a fused phase uses the smallest item
    two = [Sphere((0.3, 0.5, 0.5), 0.2), Sphere((0.75, 0.5, 0.5), 0.1)]
    nf = generate(RVE(unitcell(), Phase(:p, two; refine=InterfaceRefinement(level=1))); h=0.15).stats.nelements
    nfu = generate(RVE(unitcell(), Phase(:p, two; fuse=true, refine=InterfaceRefinement(level=1))); h=0.15).stats.nelements
    @test nfu > 1.5 * nf
    # one item fused alone is unchanged
    @test volume(geo(RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.2); fuse=true))), :p) ≈ 4π * 0.008 / 3 rtol = 1e-9
end
