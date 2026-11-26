rve = RVE(unitcell(), Phase(:p, Sphere((0.5, 0.5, 0.5), 0.2)))
bad = RVE(unitcell(), Phase(:p, Cut(Sphere((0.5, 0.5, 0.5), 0.1), Sphere((0.5, 0.5, 0.5), 0.2))))

@testset "session ownership (S1-S4)" begin
    generate(rve; geometry_only=true)
    @test gmsh.isInitialized() == 0
    @test_throws ErrorException generate(bad)
    @test gmsh.isInitialized() == 0
    # keep_open
    res = generate(rve; h=0.25, keep_open=true)
    @test gmsh.isInitialized() == 1
    f = joinpath(OUTDIR, "kept.msh")
    @test write_mesh(res, f) == f && isfile(f) && f in res.files
    close_gmsh()
    close_gmsh()
    @test gmsh.isInitialized() == 0
    # do-block
    inside = Ref(false)
    generate(rve; geometry_only=true) do r
        inside[] = gmsh.isInitialized() == 1 && r.model in gmsh.model.list()
    end
    @test inside[] && gmsh.isInitialized() == 0
    @test_throws ErrorException generate(rve; geometry_only=true) do r
        error("boom")
    end
    @test gmsh.isInitialized() == 0
    # caller-owned session stays open; a failing generate removes its model
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        generate(rve; geometry_only=true, name="first")
        @test gmsh.isInitialized() == 1 && "first" in gmsh.model.list()
        @test_throws ErrorException generate(bad; name="broken")
        @test gmsh.isInitialized() == 1 && !("broken" in gmsh.model.list())
    finally
        gmsh.finalize()
    end
end

@testset "outputs and closed results (S5)" begin
    res = generate(rve; h=0.25, output=[joinpath(OUTDIR, "a.msh"), joinpath(OUTDIR, "a.vtk")])
    @test all(isfile, res.files) && length(res.files) == 2
    @test_throws ErrorException write_mesh(res, joinpath(OUTDIR, "late.msh"))
    @test_throws ErrorException view_mesh(res)
    # view is skipped without a display
    withenv("RVES_NOPOPUP" => "1") do
        @test generate(rve; geometry_only=true, view=true) isa RVEMesh
    end
    s = sprint(show, MIME("text/plain"), res)
    @test occursin("matrix", s) && occursin("periodicity", s) && occursin("tetrahedra", s)
    @test occursin("elements", sprint(show, res))
end
