# Porous solid / open structure: no matrix, the mesh contains only the phases.
using HyperRVE

rve = RVE(Cell((1.0, 1.0, 1.0)), Phase(:solid, Sphere((0.5, 0.5, 0.5), 0.56); refine = 0.03); matrix = nothing)
res = generate(rve; h = 0.08, output = joinpath(@__DIR__, "porous.msh"))
println("solid fraction: ", volume_fraction(res, :solid), "   porosity: ", 1 - volume_fraction(res, :solid))
println(interface_names(res))     # interface_solid_void: the pore walls
