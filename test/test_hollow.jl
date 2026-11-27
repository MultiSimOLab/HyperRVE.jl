# Hollow strut lattices (tubes): shell + cavity phases.

@testset "hollow struts" begin
    c = unitcell()
    t = lattice_topology(:cubic)
    hl = hollow_struts(c, t.nodes, t.edges; radius=0.1, wall=0.03)
    @test hl.shell == struts(c, t.nodes, t.edges; radius=0.1)
    @test hl.core == struts(c, t.nodes, t.edges; radius=0.07)
    # volumes against the exact densities, hollow and solid nodes
    for node in (:hollow, :solid)
        h = hollow_lattice(c, :cubic; radius=0.1, wall=0.03, node)
        g = geo(RVE(c, [Phase(:lattice, h.shell; fuse=true), Phase(:cavity, h.core; fuse=node === :hollow, void=true)]; matrix=nothing))
        @test volume_fraction(g, :lattice) ≈ relative_density(:cubic, 0.1; wall=0.03, node) rtol = 1e-5
    end
    h = hollow_lattice(c, :bcc; radius=0.08, wall=0.03, node=:solid)
    g = geo(RVE(c, [Phase(:lattice, h.shell; fuse=true), Phase(:cavity, h.core; void=true)]; matrix=nothing))
    @test volume_fraction(g, :lattice) ≈ relative_density(:bcc, 0.08; wall=0.03, node=:solid) rtol = 1e-5
    @test relative_density(:cubic, 0.1; wall=0.03) ≈ relative_density(:cubic, 0.1) - relative_density(:cubic, 0.07)
    # wall from the density, both modes
    for node in (:hollow, :solid)
        h = hollow_lattice(c, :cubic; radius=0.1, density=0.05, node)
        w = 0.1 - first(k for k in h.core if k isa Cylinder).radius
        @test relative_density(:cubic, 0.1; wall=w, node) ≈ 0.05 rtol = 1e-8
    end
    # solid nodes: cavities stop where the struts leave the node (cubic: √2 R from the node)
    h = hollow_lattice(c, :cubic; radius=0.1, wall=0.03, node=:solid)
    @test all(k -> isapprox(sqrt(sum(abs2, k.axis)), 1 - 2 * sqrt(2) * 0.1; rtol=1e-9), h.core)
end

@testset "hollow meshes" begin
    c = unitcell()
    h = hollow_lattice(c, :cubic; radius=0.12, wall=0.04)
    res = generate(RVE(c, [Phase(:lattice, h.shell; fuse=true), Phase(:cavity, h.core; fuse=true, void=true)]; matrix=nothing); h=0.06)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :lattice) ≈ relative_density(:cubic, 0.12; wall=0.04) rtol = 1e-3
    @test res.stats.min_quality > 0
    # core-shell in a matrix: the inside is a third phase
    res = generate(RVE(c, [Phase(:shell, h.shell; fuse=true), Phase(:core, h.core; fuse=true)]); h=0.07)
    @test Set(phase_names(res)) == Set(["matrix", "shell", "core"])
    @test volume_fraction(res, :core) ≈ relative_density(:cubic, 0.08) rtol = 1e-3
end

@testset "hollow errors" begin
    c = unitcell()
    t = lattice_topology(:cubic)
    @test_throws ArgumentError hollow_struts(c, t.nodes, t.edges; radius=0.1, wall=0.1)
    @test_throws ArgumentError hollow_struts(c, t.nodes, t.edges; radius=0.1, wall=-0.01)
    @test_throws ArgumentError hollow_struts(c, t.nodes, t.edges; radius=0.1, wall=0.03, node=:open)
    @test_throws ArgumentError hollow_struts(c, t.nodes, t.edges; radius=[0.1, 0.1, 0.1], wall=0.03)
    @test_throws ArgumentError hollow_lattice(c, :cubic; radius=0.1)
    @test_throws ArgumentError hollow_lattice(c, :cubic; radius=0.1, wall=0.03, density=0.05)
    @test_throws ArgumentError hollow_lattice(c, :cubic; radius=0.1, density=0.5)
    @test_throws ArgumentError hollow_lattice(Cell((1.0, 1.0, 2.0)), :cubic; radius=0.1, density=0.05)
    @test_throws ArgumentError relative_density(:cubic, 0.1; wall=0.2)
    @test_throws ArgumentError relative_density(:square, 0.1; wall=0.02)
end

@testset "hollow review fixes" begin
    c = unitcell()
    # solid nodes on other presets and a non-cubic cell away from the origin: OCC volume = shell − cavities
    for (k, R) in ((:fccz, 0.07), (:diamond, 0.07))
        h = hollow_lattice(c, k; radius=R, wall=0.02, node=:solid)
        g = geo(RVE(c, [Phase(:lattice, h.shell; fuse=true), Phase(:cavity, h.core; void=true)]; matrix=nothing))
        @test volume_fraction(g, :lattice) ≈ relative_density(k, R; wall=0.02, node=:solid) rtol = 1e-5
    end
    cb = Cell((2.0, 1.0, 0.5); origin=(10.0, -3.0, 7.0))
    t = lattice_topology(:bcc)
    h = hollow_struts(cb, t.nodes, t.edges; radius=0.06, wall=0.02, node=:solid)
    g = geo(RVE(cb, [Phase(:lattice, h.shell; fuse=true), Phase(:cavity, h.core; void=true)]; matrix=nothing))
    solid = geo(RVE(cb, Phase(:lattice, h.shell; fuse=true); matrix=nothing))
    @test volume_fraction(g, :lattice) ≈ volume_fraction(solid, :lattice) - sum(k -> π * k.radius^2 * sqrt(sum(abs2, k.axis)), h.core) / 1.0 rtol = 1e-5
    # node radii: number and vector; hollow nodes need spheres
    tc = lattice_topology(:cubic)
    @test length(hollow_struts(c, tc.nodes, tc.edges; radius=0.1, wall=0.03, node_radius=0.15).core) > 0
    @test_throws ArgumentError hollow_struts(c, tc.nodes, tc.edges; radius=0.1, wall=0.03, node_radius=nothing)
    @test_throws ArgumentError hollow_struts(c, tc.nodes, tc.edges; radius=0.1, wall=0.03, node_radius=0.02)
    @test length(hollow_struts(c, tc.nodes, tc.edges; radius=0.1, wall=0.03, node=:solid, node_radius=nothing).core) == 3
    # short struts stay solid with a warning; none long enough is an error
    nodes = [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0)]
    edges = [(1, 2, (0, 0, 0)), (2, 1, (1, 0, 0)), (1, 1, (0, 1, 0)), (1, 1, (0, 0, 1))]
    h = @test_logs (:warn, r"too short") match_mode = :any hollow_struts(c, nodes, edges; radius=0.04, wall=0.02, node=:solid)
    @test length(h.core) == 3
    @test_throws ArgumentError hollow_lattice(c, :kelvin; radius=0.13, wall=0.03, node=:solid)
    # mesh with solid nodes: unfused cavities crossing the faces
    h = hollow_lattice(c, :bcc; radius=0.08, wall=0.03, node=:solid)
    res = generate(RVE(c, [Phase(:lattice, h.shell; fuse=true), Phase(:cavity, h.core; void=true)]; matrix=nothing); h=0.05)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :lattice) ≈ relative_density(:bcc, 0.08; wall=0.03, node=:solid) rtol = 1e-3
    # validations
    @test_throws ArgumentError relative_density(:octet, 0.1; node=:foo)
    @test_throws ArgumentError hollow_lattice(c, :cubic; radius=[0.1, 0.1, 0.1], density=0.05)
    @test_throws ArgumentError hollow_lattice(c, :square; radius=0.1, wall=0.03)
end

