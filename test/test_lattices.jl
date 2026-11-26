# Strut lattices: topology of the presets, struts(), densities, meshing.

using RVEs: lattice_topology

@testset "preset topologies (S1)" begin
    table = Dict(:cubic => (1, 3, 6), :bcc => (2, 8, 8), :fcc => (4, 12, nothing), :octet => (4, 24, 12),
                 :kelvin => (12, 24, 4), :diamond => (8, 16, 4))
    for (k, (nn, ne, val)) in table
        t = lattice_topology(k)
        @test length(t.nodes) == nn
        @test length(t.edges) == ne
        if val !== nothing
            deg = zeros(Int, nn)
            for (i, j, _) in t.edges
                deg[i] += 1; deg[j] += 1
            end
            @test all(==(val), deg)
        end
    end
    @test_throws ArgumentError lattice_topology(:nope)
end

@testset "struts() validation and emission (S2, S3)" begin
    c = unitcell()
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0, 0, 0)], [(1, 2)]; radius=0.05)
    @test_throws ArgumentError struts(c, [(0, 0, 0)], [(1, 1)]; radius=0.05)
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0.5, 0.5, 0.5)], [(1, 2), (2, 1)]; radius=0.05)
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0.5, 0.5, 0.5)], [(1, 3)]; radius=0.05)
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0.5, 0.5, 0.5)], [(1, 2)]; radius=0.6)
    @test_throws ArgumentError struts(c, [(0, 0, 0), (0.5, 0.5, 0.5)], [(1, 2)]; radius=[0.1, 0.1])
    @test_throws ArgumentError lattice(c, :bcc)
    @test_throws ArgumentError lattice(Cell((1.0, 2.0, 1.0)), :bcc; density=0.1)
    # periodic cell: struts listed once (+ node spheres)
    items = lattice(c, :bcc; radius=0.05)
    @test count(x -> x isa Cylinder, items) == 8 && count(x -> x isa Sphere, items) == 2
    # non-periodic cubic cell: the complete finite truss (12 edges, 8 corners)
    np = lattice(Cell((1.0, 1.0, 1.0); periodic=false), :cubic; radius=0.05)
    @test count(x -> x isa Cylinder, np) == 12 && count(x -> x isa Sphere, np) == 8
    # node_radius options and per-strut radii
    @test isempty(filter(x -> x isa Sphere, lattice(c, :bcc; radius=0.05, node_radius=nothing)))
    @test all(x -> x.radius == 0.07, filter(x -> x isa Sphere, lattice(c, :bcc; radius=0.05, node_radius=0.07)))
    t = lattice_topology(:bcc)
    its = struts(c, t.nodes, t.edges; radius=[0.05 + 0.005k for k in 1:8])
    @test filter(x -> x isa Sphere, its)[2].radius ≈ 0.09               # :auto = largest strut at the node
    # nodes outside [0, 1) are wrapped and offsets corrected
    @test length(struts(c, [(1.0, 1.0, 1.0), (0.5, 0.5, 0.5)], [(2, 1)]; radius=0.05)) == 3
end

@testset "densities (S4, S5)" begin
    for k in (:cubic, :bcc, :fcc, :octet, :kelvin, :diamond)
        ρ = 0.15
        x = strut_radius(k, ρ)
        @test relative_density(k, x) ≈ ρ rtol = 1e-12
    end
    @test_throws ArgumentError relative_density(:octet, 0.5)
    @test_throws ArgumentError strut_radius(:octet, 0.9)
    # exact formula vs OCC volume of the fused lattice
    for (k, x) in ((:bcc, 0.08), (:octet, 0.06), (:kelvin, 0.05), (:cubic, 0.1))
        res = geo(RVE(unitcell(), Phase(:lat, lattice(unitcell(), k; radius=x); fuse=true); matrix=nothing))
        @test volume_fraction(res, :lat) ≈ relative_density(k, x) rtol = 1e-4
    end
    # density keyword and first-order upper bound
    res = geo(RVE(unitcell(), Phase(:lat, lattice(unitcell(), :bcc; density=0.2); fuse=true); matrix=nothing))
    @test volume_fraction(res, :lat) ≈ 0.2 rtol = 1e-4
    t = lattice_topology(:bcc)
    @test relative_density(unitcell(), t.nodes, t.edges; radius=0.08) > relative_density(:bcc, 0.08)
end

@testset "meshed lattices (S6, S7)" begin
    # kept cheap: one preset with the matrix, two without (the others are checked as geometry above)
    for (k, x, m) in ((:bcc, 0.08, :matrix), (:octet, 0.06, nothing), (:diamond, 0.06, nothing))
        res = generate(RVE(unitcell(), Phase(:lat, lattice(unitcell(), k; radius=x); fuse=true); matrix=m); h=0.1)
        @test all(c -> c.ok && c.npairs > 0, res.periodicity)
        @test length(res[:lat].entities) >= 1 && res.stats.nelements > 0
    end
    # finite non-periodic block, 2 x 1 x 1
    blk = generate(RVE(Cell((1.0, 1.0, 1.0); periodic=false), Phase(:lat, lattice(Cell((1.0, 1.0, 1.0); periodic=false), :bcc; radius=0.08); fuse=true);
                       matrix=nothing, tiling=(2, 1, 1)); h=0.1)
    @test length(blk[:lat].entities) == 1 && isempty(blk.periodicity)
end
