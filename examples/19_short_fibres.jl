# Random short fibres: capsules of finite length with a prescribed orientation, periodic and
# without overlaps. `length` is the distance between the centres of the two caps.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))

# isotropic orientations
f = random_capsules(cell; length = 0.3, radius = 0.04, fraction = 0.025, rng = 1)
res = generate(RVE(cell, Phase(:fibres, f; refine = 0.03)); h = 0.1, output = joinpath(@__DIR__, "short_fibres.msh"))
A = orientation_tensor(f)
println(length(f), " isotropic fibres, fraction ", round(volume_fraction(res, :fibres); digits = 4),
        ", orientation tensor diagonal ", round.((A[1, 1], A[2, 2], A[3, 3]); digits = 2))

# random in the xy plane, in a flat cell
flat = Cell((1.0, 1.0, 0.4))
f = random_capsules(flat; length = 0.35, radius = 0.04, fraction = 0.05, orientation = :planar, axis = 3, rng = 2)
res = generate(RVE(flat, Phase(:fibres, f; refine = 0.03)); h = 0.08, output = joinpath(@__DIR__, "short_fibres_planar.msh"))
println(length(f), " in-plane fibres, fraction ", round(volume_fraction(res, :fibres); digits = 4))

# aligned with x: higher fractions are reachable
f = random_capsules(cell; length = 0.4, radius = 0.05, fraction = 0.08, orientation = (1, 0, 0), rng = 3)
res = generate(RVE(cell, Phase(:fibres, f; refine = 0.03)); h = 0.1, output = joinpath(@__DIR__, "short_fibres_aligned.msh"))
println(length(f), " aligned fibres, fraction ", round(volume_fraction(res, :fibres); digits = 4))

# hybrid composite: particles first, fibres around them
p = random_spheres(cell; radius = 0.12, n = 6, rng = 4)
f = random_capsules(cell; length = 0.3, radius = 0.04, n = 12, avoid = p, min_gap = 0.02, rng = 5)
rve = RVE(cell, [Phase(:particles, p; refine = 0.04), Phase(:fibres, f; refine = 0.03)]; overlaps = :error)
res = generate(rve; h = 0.1, output = joinpath(@__DIR__, "short_fibres_hybrid.msh"))
println(length(p), " particles + ", length(f), " fibres, fractions ", volume_fractions(res))
