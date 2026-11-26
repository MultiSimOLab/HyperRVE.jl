faces = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
edges = vcat([["edge_$(a)_$(b)" for a in ("$(p)min", "$(p)max") for b in ("$(q)min", "$(q)max")]
              for (p, q) in (("y", "z"), ("x", "z"), ("x", "y"))]...)
corners = ["corner_x$(i)_y$(j)_z$(k)" for k in ("min", "max") for j in ("min", "max") for i in ("min", "max")]

@testset "names, tags and contents (P1, P2)" begin
    cell = Cell((2.0, 1.0, 0.5); origin=(1.0, 1.0, 1.0))
    rve = RVE(cell, Phase(:p, Sphere(cellpoint(cell, (0, 0, 0)), 0.2)))
    generate(rve; geometry_only=true) do res
        @test sort(group_names(res; dim=3)) == ["matrix", "p"]
        @test group_tag(res, :matrix; dim=3) == 1 && group_tag(res, :p; dim=3) == 2
        cellgroups(dim) = Set(filter(!startswith("interface"), group_names(res; dim=dim)))
        @test cellgroups(2) == Set(vcat(faces, "boundary"))
        @test cellgroups(1) == Set(vcat(faces, "boundary", edges))
        @test cellgroups(0) == Set(vcat(faces, "boundary", edges, "corners", corners))
        for (i, f) in enumerate(faces), dim in 0:2
            @test group_tag(res, f; dim=dim) == i
        end
        @test group_tag(res, "corners"; dim=0) == 7
        @test group_tag(res, "boundary"; dim=2) == 8
        @test [group_tag(res, e; dim=1) for e in edges] == collect(11:22)
        @test [group_tag(res, c; dim=0) for c in corners] == collect(31:38)
        g(name, dim) = only(x for x in res.groups if x.name == name && x.dim == dim)
        @test length(g("corners", 0).entities) == 8
        for c in corners
            @test length(g(c, 0).entities) == 1
        end
        # the CAD model is built in a normalised frame (origin at the cell corner, longest
        # edge = 1): map back with the frame
        fr = RVEs.Frame(rve)
        s = fr.scale
        # the corner xmax_ymin_zmin is at (3, 1, 1)
        p = g("corner_xmax_ymin_zmin", 0).entities[1]
        @test all(fr.origin .+ s .* Tuple(gmsh.model.getValue(0, p, Float64[])) .≈ (3.0, 1.0, 1.0))
        # face areas, edge lengths and the whole boundary
        area(name) = s^2 * sum(gmsh.model.occ.getMass(2, t) for t in g(name, 2).entities)
        @test area("xmin") ≈ 0.5 rtol = 1e-9
        @test area("ymax") ≈ 1.0 rtol = 1e-9
        @test area("zmin") ≈ 2.0 rtol = 1e-9
        @test area("boundary") ≈ 2 * (0.5 + 1.0 + 2.0) rtol = 1e-9
        len(name) = s * sum(gmsh.model.occ.getMass(1, t) for t in g(name, 1).entities)
        @test len("edge_ymin_zmin") ≈ 2.0 rtol = 1e-9
        @test len("edge_xmax_zmax") ≈ 1.0 rtol = 1e-9
        @test len("edge_xmin_ymax") ≈ 0.5 rtol = 1e-9
        @test_throws KeyError group_tag(res, :nope; dim=3)
    end
end

@testset "file round trip (P3)" begin
    f = joinpath(OUTDIR, "roundtrip.msh")
    rve = RVE(unitcell(), [Phase(:a, Sphere((0.5, 0.5, 0.5), 0.2)), Phase(:b, Fibre((0, 0, 0), (0, 0, 1), 0.1))])
    res = generate(rve; h=0.2, output=f)
    @test res.files == [f] && isfile(f)
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(f)
        names3 = Set(gmsh.model.getPhysicalName(3, t) for (_, t) in gmsh.model.getPhysicalGroups(3))
        @test names3 == Set(["matrix", "a", "b"])
        @test "corners" in Set(gmsh.model.getPhysicalName(0, t) for (_, t) in gmsh.model.getPhysicalGroups(0))
        nel = 0
        for (_, t) in gmsh.model.getPhysicalGroups(3), e in gmsh.model.getEntitiesForPhysicalGroup(3, t)
            _, tags, _ = gmsh.model.mesh.getElements(3, e)
            nel += sum(length, tags; init=0)
        end
        @test nel == res.stats.nelements
    finally
        gmsh.finalize()
    end
end
