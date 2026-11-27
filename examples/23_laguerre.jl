# Dense polydisperse packings and Laguerre (power) tessellations: grains with a size
# distribution, as a polycrystal and as the cells of a closed-cell foam.
using HyperRVE, Random

cell = Cell((1.0, 1.0, 1.0))
# radii from a truncated log-normal law (diameters 0.25-0.6, mean about 0.35)
lognormal(rng) = (while true
    d = 0.35 * exp(0.33 * randn(rng) - 0.05); 0.25 <= d <= 0.6 && return d / 2
end)
s = dense_spheres(cell; radius = lognormal, fraction = 0.55, rng = 1)     # far denser than random_spheres
t = laguerre(cell; spheres = s, volumes = :spheres)                        # grain volumes ∝ sphere volumes
d = grain_diameters(t)
println(length(s), " grains, packing ", round(volume_fraction(s, cell); digits = 3),
        ", grain diameters ", round(minimum(d); digits = 3), "-", round(maximum(d); digits = 3))
res = generate(RVE(cell, grain_phases(t); matrix = nothing); h = 0.12, mesh = MeshOptions(optimize = :netgen),
               output = joinpath(@__DIR__, "laguerre_polycrystal.msh"))
println(length(res.phases), " grain phases, periodic: ", all(x -> x.ok, res.periodicity))

# closed-cell foam on the Laguerre cells (relative density 0.15)
w = wall_thickness(t, 0.15)
res = generate(RVE(cell, Phase(:pores, pores(t; wall = w); void = true); matrix = :walls); h = 0.08,
               mesh = MeshOptions(optimize = :netgen),
               output = joinpath(@__DIR__, "laguerre_foam.msh"))
println("wall ", round(w; digits = 4), ", solid fraction ", round(volume_fraction(res, :walls); digits = 3))

# 2D: disks packed to 0.75 and their Laguerre cells
c2 = Cell((1.0, 1.0))
t2 = laguerre(c2; spheres = dense_spheres(c2; radius = r -> 0.05 + 0.05rand(r), fraction = 0.75, rng = 2), volumes = :spheres)
res = generate(RVE(c2, grain_phases(t2); matrix = nothing); h = 0.03, output = joinpath(@__DIR__, "laguerre_2d.msh"))
println("2D: ", length(res.phases), " grains")
