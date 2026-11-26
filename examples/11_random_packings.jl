# Random packings: particles at a target volume fraction, a size distribution, fibres and
# disks. A seed (rng = ...) makes each packing reproducible.
using RVEs

cell = Cell((1.0, 1.0, 1.0))

# equal spheres, 20 % by volume, with the default gaps to each other and to the cell faces
s = random_spheres(cell; radius = 0.08, fraction = 0.2, rng = 1)
res = generate(RVE(cell, Phase(:particles, s; refine = 0.02)); h = 0.07, output = joinpath(@__DIR__, "random_spheres.msh"))
println(length(s), " spheres, packing ", round(volume_fraction(s, cell); digits = 4),
        ", meshed ", round(volume_fraction(res, :particles); digits = 4))

# lognormal radii, 15 %, kept inside the cell (no cut particles)
s = random_spheres(cell; radius = rng -> exp(log(0.06) + 0.25 * randn(rng)), fraction = 0.15,
                   boundary = :inside, rng = 2)
println(length(s), " polydisperse spheres, radii ", round(s[end].radius; digits = 3), "-", round(s[1].radius; digits = 3))

# aligned fibres along x, 40 %, in a long cell (dense packing: h must resolve the thin gaps)
fc = Cell((2.0, 0.5, 0.5))
f = random_fibres(fc; axis = 1, radius = 0.05, fraction = 0.4, rng = 3)
res = generate(RVE(fc, Phase(:fibres, f)); h = 0.03, output = joinpath(@__DIR__, "random_fibres.msh"))
println(length(f), " fibres, meshed fraction ", round(volume_fraction(res, :fibres); digits = 4))

# 2D: disks of two sizes
c2 = Cell((1.0, 1.0))
d = random_disks(c2; radius = [fill(0.08, 10); fill(0.04, 30)], rng = 4)
res = generate(RVE(c2, Phase(:disks, d)); h = 0.03, output = joinpath(@__DIR__, "random_disks.msh"))
println(length(d), " disks, fraction ", round(volume_fraction(res, :disks); digits = 4))
