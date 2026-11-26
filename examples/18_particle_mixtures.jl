# Mixtures of particle populations and clustered particles. Each population is packed around
# the ones already placed (avoid = ...) and goes to its own phase.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))

# hybrid magnetorheological elastomer: large hard-magnetic and small soft-magnetic particles.
# Pack the large population first; fractions refer to the whole cell.
hard = random_spheres(cell; radius = 0.14, fraction = 0.10, rng = 1)
soft = random_spheres(cell; radius = 0.09, fraction = 0.06, avoid = hard, rng = 2)
rve = RVE(cell, [Phase(:hard, hard), Phase(:soft, soft)]; overlaps = :error)
res = generate(rve; h = 0.1, output = joinpath(@__DIR__, "mixture_hard_soft.msh"))
println(length(hard), " hard + ", length(soft), " soft particles, meshed fractions ",
        round(volume_fraction(res, :hard); digits = 4), " and ", round(volume_fraction(res, :soft); digits = 4))

# chains of hard particles (field-structured) with soft particles in between
chains = particle_chains(cell; radius = 0.09, gap = 0.02, n_chains = 2, rng = 3)
filler = random_spheres(cell; radius = 0.06, fraction = 0.04, avoid = chains, rng = 4)
res = generate(RVE(cell, [Phase(:chains, chains), Phase(:filler, filler)]; overlaps = :error);
               h = 0.1, output = joinpath(@__DIR__, "mixture_chains_filler.msh"))
println(length(chains), " chained + ", length(filler), " free particles")

# 2D: clustered disks (5 clusters of radius 0.3) and a second, uniform population around them
c2 = Cell((2.0, 2.0))
clustered = random_disks(c2; radius = 0.04, fraction = 0.15, cluster = (n = 5, radius = 0.3), rng = 5)
uniform = random_disks(c2; radius = 0.025, fraction = 0.05, avoid = clustered, rng = 6)
res = generate(RVE(c2, [Phase(:clustered, clustered), Phase(:uniform, uniform)]; overlaps = :error); h = 0.04,
               output = joinpath(@__DIR__, "mixture_clusters_2d.msh"))
println(length(clustered), " clustered + ", length(uniform), " uniform disks, fractions ",
        round(volume_fraction(res, :clustered); digits = 4), " and ", round(volume_fraction(res, :uniform); digits = 4))
