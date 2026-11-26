# Periodic Voronoi polycrystals: one phase per grain (grain boundaries as interface groups),
# grain groups as two phases, and a 2D tessellation.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))
t = voronoi(cell; n = 20, lloyd = 2, rng = 1)          # 20 grains, two Lloyd iterations
println(t, ", grain volumes ", round(minimum(grain_volumes(t)); digits = 4), "-", round(maximum(grain_volumes(t)); digits = 4))
res = generate(RVE(cell, grain_phases(t); matrix = nothing); h = 0.08, output = joinpath(@__DIR__, "polycrystal.msh"))
println(length(res.phases), " grain phases, ", length(interface_names(res)), " grain-boundary groups")

# two-phase microstructure: grains assigned at random to phase 1 or 2
groups = [isodd(i) ? 1 : 2 for i in 1:20]
res = generate(RVE(cell, grain_phases(t; prefix = "phase", groups); matrix = nothing); h = 0.08,
               output = joinpath(@__DIR__, "two_phase_grains.msh"))
println("phase fractions: ", [round(p.fraction; digits = 3) for p in res.phases])

# 2D
c2 = Cell((1.0, 1.0))
t2 = voronoi(c2; n = 30, rng = 2)
res = generate(RVE(c2, grain_phases(t2); matrix = nothing); h = 0.03, output = joinpath(@__DIR__, "polycrystal_2d.msh"))
println("2D: ", length(res.phases), " grains")
