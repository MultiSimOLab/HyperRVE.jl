# Coated particle: phases painted in order (the later phase wins where they overlap), and
# interface physical groups for cohesive or jump conditions.
using HyperRVE

c = (0.5, 0.5, 0.5)
rve = RVE(Cell((1.0, 1.0, 1.0)), [Phase(:coating, Sphere(c, 0.3)), Phase(:core, Sphere(c, 0.25))])
res = generate(rve; h = 0.08, output = joinpath(@__DIR__, "coated.msh"))
println(interface_names(res))
println("core/coating interface: area ", interface_area(res, :coating, :core),
        ", physical tag ", interface_tag(res, :coating, :core))
# in HyperFEM / Gridap:
#   model = GmshDiscreteModel("coated.msh")
#   Γ = BoundaryTriangulation(model; tags = "interface_coating_core")
