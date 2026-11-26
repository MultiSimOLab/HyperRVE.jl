# Capsule primitive and fused phases.

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
        res = geo(RVE(unitcell(), Phase(:c, Capsule(p, q, 0.1))))
        @test volume(res, :c) ≈ capsule_vol(0.1, L) rtol = 1e-9
    end
    # crossing a face and a corner: completed by its images
    res = generate(RVE(unitcell(), Phase(:c, Capsule((0.8, 0.9, 0.5), (1.2, 1.1, 0.5), 0.08))); h=0.1)
    @test volume(res, :c) ≈ capsule_vol(0.08, sqrt(0.16 + 0.04)) rtol = 1e-6
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    # scaled, shifted cell
    k = 1e-3
    ck = Cell((k, k, k); origin=(2k, 0, 0))
    rk = geo(RVE(ck, Phase(:c, Capsule(cellpoint(ck, (0.3, 0.5, 0.5)), cellpoint(ck, (0.7, 0.5, 0.5)), 0.1k))))
    @test volume_fraction(rk, :c) ≈ capsule_vol(0.1, 0.4) rtol = 1e-9
end

@testset "fused phases (K2-K5)" begin
    r = 0.08
    un = geo(RVE(unitcell(), Phase(:lat, bcc_items(r)); matrix=nothing))
    fu = geo(RVE(unitcell(), Phase(:lat, bcc_items(r); fuse=true); matrix=nothing))
    @test volume(fu, :lat) ≈ volume(un, :lat) rtol = 1e-6
    @test length(fu[:lat].entities) == 1                          # one connected body
    @test length(un[:lat].entities) > 10
    @test occursin("fused", sprint(show, Phase(:lat, bcc_items(r); fuse=true)))
    # no internal surfaces: every surface of the open fused lattice has one adjacent volume
    generate(RVE(unitcell(), Phase(:lat, bcc_items(r); fuse=true); matrix=nothing); geometry_only=true) do res
        for (_, s) in gmsh.model.getEntities(2)
            @test length(gmsh.model.getAdjacencies(2, s)[1]) == 1
        end
    end
    # meshed, periodic, with and without matrix
    for m in (:matrix, nothing)
        res = generate(RVE(unitcell(), Phase(:lat, bcc_items(r); fuse=true); matrix=m); h=0.08)
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
        @test res.stats.min_quality > 0.05
    end
    # supercell: unit path (fuse one cell, clip, tile) equals the volume of n copies
    s2 = geo(supercell(RVE(unitcell(), Phase(:lat, bcc_items(r); fuse=true); matrix=nothing), 2))
    @test volume(s2, :lat) ≈ 8 * volume(fu, :lat) rtol = 1e-6
    @test length(s2[:lat].entities) == 1
    # non-periodic tiling: flat path; duplicates of shared boundary struts are dropped
    np = Cell((1.0, 1.0, 1.0); periodic=(false, true, true))
    f2 = geo(RVE(np, Phase(:lat, bcc_items(r); fuse=true); matrix=nothing, tiling=(2, 1, 1)))
    @test length(f2[:lat].entities) == 1
    # interfaces of a fused phase in a matrix: only matrix|lat, no same-phase surfaces
    gi = geo(RVE(unitcell(), Phase(:lat, bcc_items(r); fuse=true)))
    @test interface_names(gi) == ["interface_matrix_lat"]
end

@testset "fuse in 2D and with refinement (K6-K8)" begin
    disks = [Disk((0.4, 0.5), 0.15), Disk((0.6, 0.5), 0.15), Disk((0.0, 0.0), 0.2)]
    un = geo(RVE(Cell((1.0, 1.0)), Phase(:p, disks)))
    fu = geo(RVE(Cell((1.0, 1.0)), Phase(:p, disks; fuse=true)))
    @test volume(fu, :p) ≈ volume(un, :p) rtol = 1e-9
    @test length(fu[:p].entities) < length(un[:p].entities)
    res = generate(RVE(unitcell(), Phase(:lat, bcc_items(0.08); fuse=true, refine=InterfaceRefinement(level=2))); h=0.1)
    @test all(x -> x.ok, res.periodicity) && res.stats.nelements > 0
    # one item fused alone is unchanged
    @test volume(geo(RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.2); fuse=true))), :p) ≈ 4π * 0.008 / 3 rtol = 1e-9
end
