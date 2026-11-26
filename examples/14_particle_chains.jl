# Field-structured magnetorheological elastomer: particles aligned in chains along z, with the
# isotropic reference at the same volume fraction.
using RVEs

cell = Cell((1.0, 1.0, 1.0))
r, g = 0.08, 0.02                                      # particle radius, gap inside a chain (0.25 r)

chains = particle_chains(cell; axis = 3, radius = r, gap = g, fraction = 0.12, rng = 1)
iso = random_spheres(cell; radius = r, n = length(chains), rng = 1)
println(length(chains), " particles, fraction ", round(volume_fraction(chains, cell); digits = 4))

# the ligaments between particles of a chain need about two elements across the gap
res = generate(RVE(cell, Phase(:particles, chains; refine = 0.6g)); h = 0.1,
               output = joinpath(@__DIR__, "mre_chains.msh"))
println("chains: ", res.stats.nelements, " elements, min quality ", round(res.stats.min_quality; digits = 3))
res = generate(RVE(cell, Phase(:particles, iso; refine = 0.02)); h = 0.1, output = joinpath(@__DIR__, "mre_isotropic.msh"))

# variants: wavy and broken chains, a body-centred tetragonal arrangement, elongated particles
wavy = particle_chains(cell; radius = r, gap = g, fraction = 0.1, waviness = 0.3, rng = 2)
broken = particle_chains(cell; radius = r, gap = g, fraction = 0.15, occupancy = 0.8, rng = 3)
bct = particle_chains(cell; radius = r, gap = g, n_chains = (2, 2), arrangement = :square, stagger = :alternate, rng = 4)
needles = particle_chains(cell; radius = 0.04, aspect = 2.0, gap = g, fraction = 0.1, rng = 5)
println(length(wavy), " / ", length(broken), " / ", length(bct), " / ", length(needles), " particles")

# 2D (plane strain): chains of disks along y
c2 = Cell((1.0, 1.0))
d = particle_chains(c2; axis = 2, radius = 0.05, gap = 0.01, fraction = 0.25, rng = 6)
generate(RVE(c2, Phase(:particles, d; refine = 0.005)); h = 0.04, output = joinpath(@__DIR__, "mre_chains_2d.msh"))
