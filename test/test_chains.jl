# Particle chains (field-structured MREs).

center(x) = x isa Union{Sphere,Disk} ? x.center : x.center
rad(x) = x isa Union{Sphere,Disk} ? x.radius : minimum(x.radii)

"Minimum surface distance between particles (minimum image along periodic directions)."
function min_surface_gap(items, cell)
    g = Inf
    for a in eachindex(items), b in a+1:length(items)
        δ = collect(center(items[a]) .- center(items[b]))
        for d in eachindex(δ)
            cell.periodic[d] && (δ[d] -= cell.size[d] * round(δ[d] / cell.size[d]))
        end
        g = min(g, sqrt(sum(abs2, δ)) - rad(items[a]) - rad(items[b]))
    end
    return g
end

"Particles grouped by their transverse position (straight chains)."
function chains_of(items, axis)
    groups = Dict{Any,Vector{Any}}()
    for x in items
        k = round.(Tuple(center(x)[d] for d in eachindex(center(x)) if d != axis); digits=9)
        push!(get!(groups, k, Any[]), x)
    end
    return collect(values(groups))
end

@testset "straight chains" begin
    c = unitcell()
    r, g = 0.06, 0.012
    items = particle_chains(c; radius=r, gap=g, fraction=0.2, rng=1)
    @test all(x -> x isa Sphere && x.radius == r, items)
    @test min_surface_gap(items, c) >= g - 1e-12
    ch = chains_of(items, 3)
    m = floor(Int, 1 / (2r + g))
    @test all(k -> length(k) == m, ch)
    p = 1 / m                                         # fitted pitch: the chain closes periodically
    for k in ch
        z = sort([center(x)[3] for x in k])
        @test all(diff(z) .≈ p) && (z[1] + 1 - z[end]) ≈ p
    end
    @test volume_fraction(items, c) ≈ length(ch) * m * 4π * r^3 / 3
    @test abs(volume_fraction(items, c) - 0.2) <= m * 4π * r^3 / 3
    @test particle_chains(c; radius=r, gap=g, n_chains=5, rng=3) == particle_chains(c; radius=r, gap=g, n_chains=5, rng=3)
    @test length(particle_chains(c; radius=r, gap=g, n_chains=5, rng=3)) == 5m
    # axis and fit
    items = particle_chains(Cell((0.5, 2.0, 0.5)); axis=2, radius=0.05, gap=0.01, n_chains=3, rng=2)
    @test length(unique(round(center(x)[1]; digits=9) for x in items)) == 3
    @test_logs (:warn, r"gap grows") particle_chains(Cell((1.0, 1.0, 0.33)); radius=0.06, gap=0.012, n_chains=2, rng=1)
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0.012, n_chains=2, fit=:strict)
    @test length(particle_chains(Cell((1.0, 1.0, 0.36)); radius=0.06, gap=0.06, n_chains=2, fit=:strict, face_gap=0.01, rng=1)) == 4
    @test_throws ArgumentError particle_chains(Cell((1.0, 1.0, 0.36)); radius=0.06, gap=0.06, n_chains=2, rng=1)   # no axial offset fits the bands
    # free axis: whole chains inside, exact gap
    cf = Cell((1.0, 1.0, 1.0); periodic=(true, true, false))
    items = particle_chains(cf; radius=r, gap=g, n_chains=4, rng=1)
    @test all(x -> r <= center(x)[3] <= 1 - r, items)
    for k in chains_of(items, 3)
        @test all(diff(sort([center(x)[3] for x in k])) .≈ 2r + g)
    end
end

@testset "regular arrangements, stagger, occupancy, waviness, aspect, 2D" begin
    c = unitcell()
    sq = particle_chains(c; radius=0.06, gap=0.012, n_chains=(3, 3), arrangement=:square, stagger=0.5)
    xy = sort(unique([round.(center(x)[1:2]; digits=9) for x in sq]))
    @test xy == sort([round.(((i - 0.5) / 3, (j - 0.5) / 3); digits=9) for i in 1:3 for j in 1:3])
    alt = particle_chains(c; radius=0.06, gap=0.012, n_chains=(2, 2), arrangement=:square, stagger=:alternate, rng=1)
    z0 = Dict(round.(center(x)[1:2]; digits=9) => minimum(center(y)[3] for y in alt if round.(center(y)[1:2]; digits=9) == round.(center(x)[1:2]; digits=9)) for x in alt)
    p = 1 / floor(Int, 1 / 0.132)
    @test abs(mod(z0[(0.25, 0.25)] - z0[(0.75, 0.25)], p) - p / 2) < 1e-9
    @test_logs (:warn, r"hexagonal") particle_chains(c; radius=0.05, gap=0.01, n_chains=(3, 3), arrangement=:hex)
    @test_throws ArgumentError particle_chains(c; radius=0.2, gap=0.01, n_chains=(3, 3), arrangement=:square)
    # occupancy
    full = particle_chains(c; radius=0.03, gap=0.006, n_chains=60, rng=4)
    part = particle_chains(c; radius=0.03, gap=0.006, n_chains=60, occupancy=0.7, rng=4)
    @test abs(length(part) / length(full) - 0.7) < 4 * sqrt(0.21 / length(full))
    # waviness: gaps never shrink, chains stay within ±w r of their axis
    w = particle_chains(c; radius=0.06, gap=0.012, fraction=0.12, waviness=0.3, rng=1)
    @test min_surface_gap(w, c) >= 0.012 - 1e-12
    # spheroids along the chain
    e = particle_chains(c; radius=0.04, aspect=2.0, gap=0.012, n_chains=3, rng=1)
    @test all(x -> x isa Ellipsoid && x.radii == (0.04, 0.04, 0.08), e)
    @test volume_fraction(e, c) ≈ length(e) * 4π * 0.04^2 * 0.08 / 3
    # 2D
    c2 = Cell((1.0, 1.0))
    d = particle_chains(c2; axis=2, radius=0.05, gap=0.005, fraction=0.25, rng=3)
    @test all(x -> x isa Disk, d) && min_surface_gap(d, c2) >= 0.005 - 1e-12
    @test length(particle_chains(c2; radius=0.05, gap=0.01, n_chains=4, arrangement=:regular)) == 4 * floor(Int, 1 / 0.11)
    # errors
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0, n_chains=2)
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0.01)
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0.01, n_chains=2, fraction=0.1)
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0.01, n_chains=2, axis=4)
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0.01, n_chains=2, occupancy=0)
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0.01, n_chains=2, arrangement=:grid)
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0.01, n_chains=2, stagger=:alternate)
    @test_throws ArgumentError particle_chains(c; radius=0.06, gap=0.01, n_chains=2, waviness=0.2, aspect=2.0)
    @test_throws ErrorException particle_chains(c; radius=0.06, gap=0.012, fraction=0.45, rng=1, max_tries=2000)
end

@testset "meshed chains" begin
    c = unitcell()
    items = particle_chains(c; radius=0.1, gap=0.03, n_chains=2, rng=2)
    res = generate(RVE(c, Phase(:particles, items; refine=0.015)); h=0.1)
    @test all(x -> x.ok && x.npairs > 0, res.periodicity)
    @test volume(res, :particles) ≈ volume_fraction(items, c) rtol = 1e-6
    @test "interface_matrix_particles" in interface_names(res)
    c2 = Cell((1.0, 1.0))
    res = generate(RVE(c2, Phase(:p, particle_chains(c2; radius=0.06, gap=0.012, n_chains=3, rng=1); refine=0.006)); h=0.05)
    @test all(x -> x.ok, res.periodicity)
end
