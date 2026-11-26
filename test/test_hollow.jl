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
        w = 0.1 - (node === :hollow ? h.core[1].radius : h.core[1].radius)
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
