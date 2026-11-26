# Inclusion built with boolean operations: a sphere clipped by a box, with a cross of fibres removed.
using RVEs

cell = Cell((1.0, 1.0, 1.0))
c = cellpoint(cell, (0.5, 0.5, 0.5))
L = 0.38
cross = Fuse(Fibre(c, (1, 0, 0), 0.1), Fibre(c, (0, 1, 0), 0.1), Fibre(c, (0, 0, 1), 0.1))
core = Intersect(Sphere(c, 0.25), Box(center = c, size = (L, L, L)))
rve = RVE(cell, Phase(:inclusion, Cut(core, cross); refine = 0.02))

res = generate(rve; h = 0.1, output = joinpath(@__DIR__, "booleans.msh"))
display(res)
