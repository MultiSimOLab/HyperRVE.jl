# Plate lattices: closed cells by clipping, exact densities, meshes.

@testset "plate densities" begin
    for x in (0.02, 0.07, 0.15)
        @test plate_density(:sc, x) ≈ 1 - (1 - x)^3 atol = 1e-12
        @test plate_density(:bcc, x) ≈ 1 - (1 - 2sqrt(2) * x)^3 atol = 1e-12
        @test plate_density(:fcc, x) ≈ 1 - max(0, 1 - 2sqrt(3) * x)^3 / 3 - 2 / 3 * (1 - sqrt(3) * x)^3 atol = 1e-12
        @test plate_density((:sc, :bcc), (x, x / 2)) ≈ 1 - (1 - x - sqrt(2) * x)^3 atol = 1e-12
    end
    # first-order estimates hold for thin plates
    @test plate_density(:fcc, 1e-4) ≈ 4sqrt(3) * 1e-4 rtol = 1e-3
    # thickness from density, both ways; isotropic SC + FCC ratio
    for k in (:sc, :bcc, :fcc)
        t = plate_thickness(k, 0.2)
        @test plate_density(k, t) ≈ 0.2 atol = 1e-10
    end
    t = plate_thickness((:sc, :fcc), 0.3)
    @test t[1] / t[2] ≈ 8sqrt(3) / 9
    @test plate_density((:sc, :fcc), t) ≈ 0.3 atol = 1e-10
    t = plate_thickness((:sc, :bcc), 0.25; ratio=(1, 2))
    @test t[2] ≈ 2t[1]
    @test plate_density((:sc, :bcc), t) ≈ 0.25 atol = 1e-10
    # any box cell
    box = Cell((1.0, 0.8, 1.5); origin=(0.2, 0.0, -0.3))
    t = plate_thickness(box, :fcc, 0.15)
    @test plate_density(box, :fcc; thickness=t) ≈ 0.15 atol = 1e-10
end

@testset "plate cells" begin
    c = Cell((1.0, 1.0, 1.0))
    for (k, ncells) in ((:sc, 1), (:bcc, 12), (:fcc, 12), ((:sc, :bcc), 24), ((:sc, :fcc), 21))
        t = k isa Symbol ? 0.05 : (0.05, 0.03)
        cells = plate_cells(c, k; thickness=t)
        @test length(cells) == ncells && all(x -> x isa Polyhedron, cells)
        @test 1 - sum(x -> abs(HyperRVE._polyhedron_volume(x.vertices, x.faces)), cells) ≈ plate_density(c, k; thickness=t) atol = 1e-12
        # OCC agrees: the matrix (plates) has the exact density
        g = generate(RVE(c, Phase(:cells, cells; void=true); matrix=:plates); geometry_only=true)
        @test volume_fraction(g, :plates) ≈ plate_density(c, k; thickness=t) rtol = 1e-7
    end
    # by density; non-cubic cell away from the origin; free faces; tiling
    box = Cell((1.0, 0.8, 1.5); origin=(0.2, 0.0, -0.3))
    cells = plate_cells(box, (:sc, :fcc); density=0.3)
    g = generate(RVE(box, Phase(:cells, cells; void=true); matrix=:plates); geometry_only=true)
    @test volume_fraction(g, :plates) ≈ 0.3 rtol = 1e-6
    cf = Cell((1.0, 1.0, 1.0); periodic=(true, true, false))
    cells = plate_cells(cf, :sc; thickness=0.1)
    g = generate(RVE(cf, Phase(:cells, cells; void=true); matrix=:plates); geometry_only=true)
    @test volume_fraction(g, :plates) ≈ 1 - 0.9^2 * 0.9 rtol = 1e-7       # half plates on the free faces
    cells = plate_cells(c, :fcc; thickness=0.04)
    g = generate(RVE(c, Phase(:cells, cells; void=true); matrix=:plates, tiling=2); geometry_only=true)
    @test volume_fraction(g, :plates) ≈ plate_density(:fcc, 0.04) rtol = 1e-7
end

@testset "plate meshes" begin
    c = Cell((1.0, 1.0, 1.0))
    for k in (:fcc, (:sc, :fcc))
        cells = plate_cells(c, k; density=0.25)
        res = generate(RVE(c, Phase(:cells, cells; void=true); matrix=:plates); h=0.06)
        @test all(x -> x.ok && x.npairs > 0, res.periodicity)
        @test volume_fraction(res, :plates) ≈ 0.25 rtol = 1e-3
        @test res.stats.min_quality > 0
    end
    # cells filled with a second material: plates and filler with their interface
    cells = plate_cells(c, :sc; thickness=0.1)
    res = generate(RVE(c, Phase(:filler, cells); matrix=:plates); h=0.08)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :plates) ≈ 1 - 0.9^3 rtol = 1e-3
    @test "interface_plates_filler" in interface_names(res) || "interface_filler_plates" in interface_names(res)
end

@testset "plate errors" begin
    c = Cell((1.0, 1.0, 1.0))
    @test_throws ArgumentError plate_cells(c, :hcp; thickness=0.05)
    @test_throws ArgumentError plate_cells(c, (:sc, :sc); thickness=0.05)
    @test_throws ArgumentError plate_cells(c, :sc)
    @test_throws ArgumentError plate_cells(c, :sc; thickness=0.05, density=0.2)
    @test_throws ArgumentError plate_cells(c, :sc; thickness=-0.05)
    @test_throws ArgumentError plate_cells(c, (:sc, :fcc); thickness=(0.05, 0.02, 0.01))
    @test_throws ArgumentError plate_cells(c, :sc; density=1.2)
    @test_throws ArgumentError plate_cells(c, :sc; thickness=1.5)                       # the plates fill the cell
    @test_throws ArgumentError plate_thickness((:sc, :bcc), 0.2)                       # no isotropic ratio known
    @test_throws ArgumentError plate_thickness((:sc, :fcc), 0.2; ratio=(1, -1))
end

@testset "plate review fixes" begin
    # any scale: the same lattice as the unit cube
    for a in (1e-9, 1e-6, 1500.0, 1e4)
        ca = Cell((a, a, a))
        for k in (:bcc, :fcc, (:sc, :fcc))
            t = k isa Symbol ? 0.05 : (0.07, 0.04)
            @test plate_density(ca, k; thickness=t .* a) ≈ plate_density(k, t) atol = 1e-12
            cells = plate_cells(ca, k; thickness=t .* a)
            @test length(cells) == length(plate_cells(Cell((1.0, 1.0, 1.0)), k; thickness=t))
        end
    end
    @test plate_density(Cell((1.0, 1.0, 5000.0)), :fcc; thickness=0.05) > 0
    @test length(plate_cells(Cell((1.0, 1.0, 5000.0)), :fcc; thickness=0.05)) > 0
    # all families, BCC + FCC
    c = Cell((1.0, 1.0, 1.0))
    for (k, n) in (((:bcc, :fcc), 84), ((:sc, :bcc, :fcc), 96))
        t = ntuple(_ -> 0.02, length(k))
        cells = plate_cells(c, k; thickness=t)
        @test length(cells) == n
        @test 1 - sum(x -> abs(HyperRVE._polyhedron_volume(x.vertices, x.faces)), cells) ≈ plate_density(c, k; thickness=t) atol = 1e-12
    end
    # vanished FCC tetrahedra: warning, 4 cells and the closed formula
    x = 0.32
    cells = @test_logs (:warn, r"vanished") plate_cells(c, :fcc; thickness=x)
    @test length(cells) == 4
    @test plate_density(:fcc, x) ≈ 1 - 2 / 3 * (1 - sqrt(3) * x)^3 atol = 1e-12
    # critical SC + FCC ratio: clear error or cells, never a cryptic failure
    for d in (-1e-9, 1e-8, 1e-6, 1e-5)
        try
            plate_cells(c, (:sc, :fcc); thickness=(sqrt(3) * (1 + d) * 0.04, 0.04))
            @test true
        catch e
            @test e isa ArgumentError && occursin("degenerate", e.msg)
        end
    end
    # never a silently wrong lattice near the critical ratio: cells add up to the density, or a clear error
    for d in (-4.4e-9, -1e-10, 3e-9), r in ((sqrt(3) * (1 + d), 1.0), (1.7320508, 1.0), (1.73205081, 1.0))
        try
            cells = plate_cells(c, (:sc, :fcc); density=0.3, ratio=r)
            @test 1 - sum(x -> abs(HyperRVE._polyhedron_volume(x.vertices, x.faces)), cells) ≈ 0.3 atol = 1e-8
        catch e
            @test e isa ArgumentError && occursin("degenerate", e.msg)
        end
    end
    @test length(plate_cells(Cell((2000.0, 2000.0, 2000.0)), (:sc, :fcc); density=0.3)) == 21
    @test_logs (:warn, r"edges") match_mode = :any plate_cells(c, (:sc, :fcc); thickness=(sqrt(3) * 1.001 * 0.04, 0.04))
    # free directions with BCC and FCC
    cf = Cell((1.0, 1.0, 1.0); periodic=(true, true, false))
    for k in (:bcc, :fcc)
        cells = plate_cells(cf, k; thickness=0.04)
        g = generate(RVE(cf, Phase(:cells, cells; void=true); matrix=:plates); geometry_only=true)
        @test volume_fraction(g, :plates) ≈ plate_density(cf, k; thickness=0.04) rtol = 1e-7
    end
    # BCC mesh
    res = generate(RVE(c, Phase(:cells, plate_cells(c, :bcc; density=0.2); void=true); matrix=:plates); h=0.06)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume_fraction(res, :plates) ≈ 0.2 rtol = 1e-6
    # extreme densities and a non-cubic cell with a tuple and a ratio
    for ρ in (1e-4, 0.99)
        @test plate_density(:sc, plate_thickness(:sc, ρ)) ≈ ρ atol = 1e-10
    end
    box = Cell((1.0, 0.8, 1.5); origin=(0.2, 0.0, -0.3))
    t = plate_thickness(box, (:sc, :bcc), 0.2; ratio=(2, 1))
    @test t[1] ≈ 2t[2]
    @test plate_density(box, (:sc, :bcc); thickness=t) ≈ 0.2 atol = 1e-10
    t1 = plate_thickness((:sc, :fcc), 0.3); t2 = plate_thickness((:fcc, :sc), 0.3)
    @test t1[1] ≈ t2[2] && t1[2] ≈ t2[1]
    # errors
    @test_throws ArgumentError plate_cells(c, (:bcc, :fcc); density=0.2)
    @test_throws ArgumentError plate_cells(c, (:sc, :bcc, :fcc); density=0.2)
    @test_throws ArgumentError plate_cells(c, (:sc, :fcc); density=0.2, ratio=(1, 2, 3))
    @test_throws ArgumentError plate_cells(c, :sc; density=0.2, ratio=(1,))
    @test_throws ArgumentError plate_cells(c, (:sc, :fcc); thickness=(0.05, 0.02), ratio=(1, 2))
    @test_throws ArgumentError plate_cells(c, :sc; density=0.0)
    @test_throws ArgumentError plate_cells(c, :sc; density=1.0)
    @test_throws ArgumentError plate_thickness(:sc, 0.2; ratio=:foo)
    @test_throws ArgumentError plate_cells(Cell((1.0, 1.0)), :sc; thickness=0.1)
end

