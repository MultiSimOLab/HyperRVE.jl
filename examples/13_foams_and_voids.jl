# Closed-cell and open-cell foams from a Voronoi tessellation, and void phases.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))
t = voronoi(cell; n = 8, lloyd = 2, rng = 1)

# closed cells: walls with 15 % relative density, empty pores
w = wall_thickness(t, 0.15)
res = generate(RVE(cell, Phase(:pores, pores(t; wall = w); void = true); matrix = :walls); h = w,
               output = joinpath(@__DIR__, "closed_foam.msh"))
println("closed foam: wall ", round(w; digits = 4), ", density ", round(res.phases[1].fraction; digits = 4))

# open cells: the edges of a tessellation as struts (fused); BCC seeds give the Kelvin foam
kelvin = voronoi(cell; seeds = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)])
res = generate(RVE(cell, Phase(:foam, foam_struts(kelvin; radius = 0.04); fuse = true); matrix = nothing); h = 0.025,
               output = joinpath(@__DIR__, "open_foam.msh"))
println("open (Kelvin) foam: density ", round(volume_fraction(res, :foam); digits = 4))

# a void phase: spherical pores in a matrix with stiff particles
items = random_spheres(cell; radius = 0.1, n = 6, rng = 2)
res = generate(RVE(cell, [Phase(:particles, items[1:3]), Phase(:pores, items[4:6]; void = true)]); h = 0.07)
println(interface_names(res))
