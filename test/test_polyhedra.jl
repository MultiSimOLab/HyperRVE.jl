# Polyhedron and Polygon primitives: validation, traits, volumes, periodic images.

cube(c, a) = Polyhedron([c .+ a .* (x, y, z) for z in (-0.5, 0.5) for y in (-0.5, 0.5) for x in (-0.5, 0.5)],
                        [[1, 3, 4, 2], [5, 6, 8, 7], [1, 2, 6, 5], [3, 7, 8, 4], [1, 5, 7, 3], [2, 4, 8, 6]])

"Prism over a 2D polygon `poly` between heights z0 and z1."
function prism(poly, z0, z1)
    n = length(poly)
    v = [[(x, y, z0) for (x, y) in poly]; [(x, y, z1) for (x, y) in poly]]
    f = [collect(1:n), collect(2n:-1:n+1), [[k, mod1(k + 1, n), mod1(k + 1, n) + n, k + n] for k in 1:n]...]
    return Polyhedron(v, f)
end

@testset "validation and traits" begin
    p = cube((0.5, 0.5, 0.5), 0.4)
    @test all(RVEs.refpoint(p) .≈ (0.5, 0.5, 0.5))
    @test RVEs.characteristic_length(p) ≈ 0.2
    k = Polyhedron(; vertices=p.vertices, faces=p.faces)
    @test k.vertices == p.vertices && k.faces == p.faces
    # faces given with mixed orientations: same volume and centroid
    q = Polyhedron(p.vertices, [i == 2 ? reverse(f) : f for (i, f) in enumerate(p.faces)])
    @test abs(RVEs._polyhedron_volume(q.vertices, q.faces)) ≈ 0.064
    # non-convex L-shaped prism
    L = prism([(0.2, 0.2), (0.8, 0.2), (0.8, 0.4), (0.4, 0.4), (0.4, 0.8), (0.2, 0.8)], 0.3, 0.7)
    @test abs(RVEs._polyhedron_volume(L.vertices, L.faces)) ≈ (0.6 * 0.2 + 0.2 * 0.4) * 0.4
    tet = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    @test_throws ArgumentError Polyhedron(tet, [[1, 2, 3], [1, 2, 4], [2, 3, 4]])                 # too few faces
    @test_throws ArgumentError Polyhedron([tet; (1, 1, 1)], [[1, 2, 3], [1, 2, 4], [2, 3, 4], [1, 3, 5]])  # open
    @test_throws ArgumentError Polyhedron(tet, [[1, 2, 3], [1, 2, 4], [2, 3, 4], [1, 3, 9]])         # missing vertex
    @test_throws ArgumentError Polyhedron(tet, [[1, 2, 2], [1, 2, 4], [2, 3, 4], [1, 3, 4]])         # repeated vertex
    @test_throws ArgumentError Polyhedron([(0, 0, 0), (1, 0, 0), (1, 1, 0.1), (0, 1, 0), (0.5, 0.5, 1)],
                                          [[1, 2, 3, 4], [1, 2, 5], [2, 3, 5], [3, 4, 5], [4, 1, 5]])  # not planar
    @test_throws ArgumentError Polyhedron(tet, [[1, 2, 3], [1, 2, 4], [2, 3, 4], [1, 3, 4.5]])
    @test_throws ArgumentError Polygon([(0, 0), (1, 0)])
    @test_throws ArgumentError Polygon([(0, 0), (1, 0), (2, 0)])                                    # zero area
    @test_throws ArgumentError Polygon([(0, 0), (1, 0), (1, 0), (0, 1)])                             # repeated vertex
    g = Polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])
    @test RVEs.refpoint(g) == (1.0, 0.5) && RVEs.characteristic_length(g) ≈ 2 * 2 / 6
end

@testset "volumes and periodic images" begin
    cell = unitcell()
    @test volume(geo(RVE(cell, Phase(:c, cube((0.5, 0.5, 0.5), 0.4)))), :c) ≈ 0.064 rtol = 1e-10
    res = generate(RVE(cell, Phase(:c, cube((0.95, 0.9, 0.1), 0.3))); h=0.08)       # crosses a corner
    @test volume(res, :c) ≈ 0.027 rtol = 1e-10
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    hex = [(0.5 + 0.15cos(k * π / 3), 0.5 + 0.15sin(k * π / 3)) for k in 0:5]
    res = geo(RVE(cell, Phase(:h, Rotate(prism(hex, 0.2, 0.8); axis=(1, 0, 0), angle=0.3))))
    @test volume(res, :h) ≈ 6 * 0.5 * 0.15^2 * sin(π / 3) * 0.6 rtol = 1e-10
    # two touching prisms in one phase and a non-convex one in another
    L = prism([(0.2, 0.2), (0.8, 0.2), (0.8, 0.4), (0.4, 0.4), (0.4, 0.8), (0.2, 0.8)], 0.3, 0.7)
    res = geo(RVE(cell, [Phase(:L, L), Phase(:c, cube((0.6, 0.6, 0.5), 0.2))]))
    @test volume(res, :L) ≈ 0.08 rtol = 1e-10
    @test volume(res, :c) ≈ 0.008 rtol = 1e-10
    # tiny cell (normalised frame)
    ck = Cell((1e-3, 1e-3, 1e-3); origin=(0.3, 0, 0))
    @test volume_fraction(geo(RVE(ck, Phase(:c, cube((0.3005, 0.0005, 0.0005), 4e-4)))), :c) ≈ 0.064 rtol = 1e-9
    # 2D polygon across a face
    c2 = Cell((1.0, 1.0))
    res = generate(RVE(c2, Phase(:p, Polygon([(0.8, 0.1), (1.2, 0.2), (1.0, 0.5), (0.7, 0.4)]))); h=0.05)
    @test volume(res, :p) ≈ 0.12 rtol = 1e-10
    @test all(x -> x.ok, res.periodicity)
end
