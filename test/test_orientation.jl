# Prescribed second-order orientation tensors for random_capsules.

@testset "central angular Gaussian fit" begin
    @test HyperRVE._acg_moments([1, 1, 1] ./ 3) ≈ fill(1 / 3, 3) atol = 1e-9
    for λ in ([0.7, 0.25, 0.05], [0.5, 0.5, 0.0], [0.9, 0.1, 0.0], [0.4, 0.35, 0.25])
        b = HyperRVE._acg_fit(λ)
        @test HyperRVE._acg_moments(b) ≈ λ atol = 1e-8
    end
    @test HyperRVE._acg_fit([0.0, 0.0, 1.0]) == [0.0, 0.0, 1.0]
end

@testset "prescribed orientation tensor" begin
    c = Cell((3.0, 3.0, 3.0))
    tens(v) = [a == b ? v[a] : 0.0 for a in 1:3, b in 1:3]
    for A in (tens([0.7, 0.25, 0.05]), tens([0.5, 0.5, 0.0]), tens([0.0, 0.0, 1.0]))
        f = random_capsules(c; length=0.15, radius=0.004, n=1500, orientation=A, rng=1)
        @test all(abs.(orientation_tensor(f) .- A) .< 0.04)
    end
    # a rotated tensor (off-diagonal terms)
    θ = π / 5
    R = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
    A = R * [0.7 0 0; 0 0.25 0; 0 0 0.05] * R'
    f = random_capsules(c; length=0.15, radius=0.004, n=1500, orientation=A, rng=2)
    @test all(abs.(orientation_tensor(f) .- A) .< 0.04)
    # same seed, same fibres
    @test random_capsules(c; length=0.15, radius=0.004, n=20, orientation=A, rng=3) ==
          random_capsules(c; length=0.15, radius=0.004, n=20, orientation=A, rng=3)
    # errors
    u = unitcell()
    for bad in ([0.5 0.1 0; 0 0.5 0; 0 0 0], [0.5 0 0; 0 0.6 0; 0 0 0], [1.2 0 0; 0 -0.2 0; 0 0 0], [0.5 0; 0 0.5], [NaN 0 0; 0 0.5 0; 0 0 0.5])
        @test_throws ArgumentError random_capsules(u; length=0.3, radius=0.03, n=3, orientation=bad)
    end
end

@testset "orientation tensor review fixes" begin
    # nearly degenerate tensors: relative accuracy of the fit
    for λ in ([1 - 1e-6, 1e-6, 0.0], [0.5, 0.5 - 1e-10, 1e-10], [1 - 3e-13, 2e-13, 1e-13])
        A = HyperRVE._acg_moments(HyperRVE._acg_fit(λ))
        @test all(i -> abs(A[i] - λ[i]) <= 1e-8 * λ[i] + 1e-14, 1:3)
    end
    # integer, Diagonal and Symmetric input; non-numeric input and I
    c = Cell((3.0, 3.0, 3.0))
    @test length(random_capsules(c; length=0.15, radius=0.004, n=5, orientation=[1 0 0; 0 0 0; 0 0 0], rng=1)) == 5
    @test length(random_capsules(c; length=0.15, radius=0.004, n=5, orientation=HyperRVE.LinearAlgebra.Diagonal([0.5, 0.5, 0.0]), rng=1)) == 5
    @test_throws ArgumentError random_capsules(c; length=0.15, radius=0.004, n=5, orientation=fill("a", 3, 3))
    @test_throws ArgumentError random_capsules(c; length=0.15, radius=0.004, n=5, orientation=HyperRVE.LinearAlgebra.I / 3)
    # redrawn orientations are reported
    A = [1 0 0; 0 1 0; 0 0 1] ./ 3
    @test_logs (:warn, r"biased") match_mode = :any random_capsules(Cell((1.0, 1.0, 0.3)); length=0.5, radius=0.02, n=12,
                                                                    orientation=A, boundary=:inside, rng=1)
end

