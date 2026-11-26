# Interface physical groups between phases.

@testset "one sphere, crossing a face, coated particle" begin
    r = 0.2
    generate(RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), r))); h=0.2) do res
        @test interface_names(res) == ["interface_matrix_p"]
        @test interface_area(res, :matrix, :p) ≈ 4π * r^2 rtol = 1e-6
        @test interface_area(res, :p, :matrix) == interface_area(res, :matrix, :p)
        @test interface_tag(res, :matrix, :p) == 100 + 2
        @test "interface_matrix_p" in group_names(res; dim=2) && "interfaces" in group_names(res; dim=2)
        @test group_tag(res, "interfaces"; dim=2) == 99
        # a closed sphere surface has only its seam curve/points
        @test "interface_matrix_p" in group_names(res; dim=1)
    end
    generate(RVE(unitcell(), Phase(:p, Sphere((1.0, 0.5, 0.5), r))); h=0.2) do res
        @test interface_area(res, :matrix, :p) ≈ 4π * r^2 rtol = 1e-6
        g(name, dim) = only(x for x in res.groups if x.name == name && x.dim == dim)
        ifs = Set(g("interface_matrix_p", 2).entities)
        faces = Set(vcat(g("xmin", 2).entities, g("xmax", 2).entities))
        @test isempty(intersect(ifs, faces))                       # cut discs are not interface
        circles = intersect(Set(g("interface_matrix_p", 1).entities), Set(g("xmin", 1).entities))
        @test !isempty(circles)                                    # the meeting circle is in both
        nel = 0
        for e in g("interface_matrix_p", 2).entities
            _, t, _ = gmsh.model.mesh.getElements(2, e); nel += sum(length, t; init=0)
        end
        @test nel > 0
    end
    R, rc = 0.3, 0.2
    res = generate(RVE(unitcell(), [Phase(:coating, Sphere((0.5, 0.5, 0.5), R)), Phase(:core, Sphere((0.5, 0.5, 0.5), rc))]); h=0.15)
    @test Set(interface_names(res)) == Set(["interface_matrix_coating", "interface_coating_core"])
    @test interface_area(res, :matrix, :coating) ≈ 4π * R^2 rtol = 1e-6
    @test interface_area(res, :coating, :core) ≈ 4π * rc^2 rtol = 1e-6
    @test_throws KeyError interface_area(res, :matrix, :core)
    @test interface_tag(res, :matrix, :coating) == 100 + 2 && interface_tag(res, :coating, :core) == 100 + 4 + 3
end

@testset "same-phase overlaps, void, disabled groups" begin
    r, d = 0.2, 0.2                                      # two spheres, centres 0.2 apart
    res = generate(RVE(unitcell(), Phase(:p, [Sphere((0.4, 0.5, 0.5), r), Sphere((0.6, 0.5, 0.5), r)])); geometry_only=true)
    cap = 2π * r * (r - d / 2)                            # area of each removed cap
    @test interface_area(res, :matrix, :p) ≈ 2 * (4π * r^2 - cap) rtol = 1e-6
    pores = generate(RVE(unitcell(), Phase(:solid, Sphere((0.5, 0.5, 0.5), 0.3)); matrix=nothing); geometry_only=true)
    @test interface_names(pores) == ["interface_solid_void"]
    @test interface_area(pores, :solid, :void) ≈ 4π * 0.09 rtol = 1e-6
    off = generate(RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.2))); geometry_only=true, interfaces=false)
    @test !any(startswith("interface"), group_names(off))
    @test all(i -> i.tag == 0, off.interfaces) && length(off.interfaces) == 1
    @test_throws ErrorException interface_tag(off, :matrix, :p)
    @test interface_area(off, :matrix, :p) ≈ 4π * 0.04 rtol = 1e-6
end

@testset "union group, unique tags, 2D" begin
    rve = RVE(unitcell(), [Phase(:a, Sphere((0.3, 0.5, 0.5), 0.15)), Phase(:b, Sphere((0.7, 0.5, 0.5), 0.15)),
                           Phase(:c, Sphere((0.5, 0.5, 0.5), 0.1))])
    generate(rve; geometry_only=true) do res
        g(name, dim) = Set(only(x for x in res.groups if x.name == name && x.dim == dim).entities)
        u = union((g(n, 2) for n in interface_names(res))...)
        @test g("interfaces", 2) == u
        for dim in 0:3
            names = Dict{Int,String}()
            for x in res.groups
                x.dim == dim || continue
                @test get!(names, x.tag, x.name) == x.name        # one name per tag and dimension
            end
        end
    end
    r2 = generate(RVE(Cell((1.0, 1.0)), Phase(:p, Disk((0.5, 0.5), 0.2))); h=0.1)
    @test interface_area(r2, :matrix, :p) ≈ 2π * 0.2 rtol = 1e-6
    @test "interface_matrix_p" in group_names(r2; dim=1)
    @test "interface_matrix_p" in group_names(r2; dim=0)          # the seam point of the circle
    @test occursin("interfaces:", sprint(show, MIME("text/plain"), r2))
end

@testset "GridapGmsh-style reading: names in the file" begin
    f = joinpath(OUTDIR, "ifs.msh")
    generate(RVE(unitcell(), [Phase(:coating, Sphere((0.5, 0.5, 0.5), 0.3)), Phase(:core, Sphere((0.5, 0.5, 0.5), 0.2))]);
             h=0.15, output=f)
    gmsh.initialize(String[], false)
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f)
        names2 = Set(gmsh.model.getPhysicalName(d, t) for (d, t) in gmsh.model.getPhysicalGroups(2))
        @test "interface_coating_core" in names2 && "interfaces" in names2 && "xmin" in names2
    finally
        gmsh.finalize()
    end
end

@testset "interfaces: scaled cells, all dimensions, void tags, names (review)" begin
    # areas in physical units (scaled, shifted cell; 3D and 2D)
    k = 1e-3
    c3 = Cell((k, 2k, k); origin=(5k, -k, 0))
    r3 = generate(RVE(c3, Phase(:p, Sphere(cellpoint(c3, (0.5, 0.5, 0.5)), 0.3k))); geometry_only=true)
    @test interface_area(r3, :matrix, :p) ≈ 4π * (0.3k)^2 rtol = 1e-9
    c2 = Cell((2k, k); origin=(1k, 3k))
    r2 = generate(RVE(c2, Phase(:p, Disk(cellpoint(c2, (0.5, 0.5)), 0.2k))); geometry_only=true)
    @test interface_area(r2, :matrix, :p) ≈ 2π * 0.2k rtol = 1e-9
    # a sphere crossing a face: pair and union groups in dimensions 2, 1 and 0
    generate(RVE(unitcell(), Phase(:p, Sphere((1.0, 0.5, 0.5), 0.2))); geometry_only=true) do res
        for dim in 0:2
            @test "interface_matrix_p" in group_names(res; dim=dim)
            @test "interfaces" in group_names(res; dim=dim)
        end
    end
    # void: tag and union
    generate(RVE(unitcell(), [Phase(:s, Sphere((0.3, 0.5, 0.5), 0.2)), Phase(:t, Sphere((0.7, 0.5, 0.5), 0.15))]; matrix=nothing);
             geometry_only=true) do res
        @test Set(interface_names(res)) == Set(["interface_s_void", "interface_t_void"])
        @test interface_tag(res, :s, :void) == 100 + 0 * 3 + 3 && interface_tag(res, :t, :void) == 100 + 1 * 3 + 3
        g(name, dim) = Set(only(x for x in res.groups if x.name == name && x.dim == dim).entities)
        for dim in 0:2
            @test g("interfaces", dim) == union(g("interface_s_void", dim), g("interface_t_void", dim))
        end
    end
    # same-phase surfaces are not interfaces
    ov = generate(RVE(unitcell(), Phase(:p, [Sphere((0.4, 0.5, 0.5), 0.2), Sphere((0.6, 0.5, 0.5), 0.2)])); geometry_only=true)
    @test interface_names(ov) == ["interface_matrix_p"]
    # one tag per name in every dimension
    generate(RVE(unitcell(), [Phase(:a, Sphere((0.3, 0.5, 0.5), 0.15)), Phase(:b, Sphere((0.7, 0.5, 0.5), 0.15))]);
             geometry_only=true) do res
        for dim in 0:3
            tags = Dict{String,Int}()
            for x in res.groups
                x.dim == dim || continue
                @test get!(tags, x.name, x.tag) == x.tag
            end
        end
    end
    # ambiguous interface names are rejected
    @test_throws ArgumentError RVE(unitcell(), [Phase(:b_c, Sphere((0.2, 0.5, 0.5), 0.1)), Phase(:a_b, Sphere((0.5, 0.5, 0.5), 0.1)),
                                              Phase(:c, Sphere((0.8, 0.5, 0.5), 0.1))]; matrix=:a)
end
