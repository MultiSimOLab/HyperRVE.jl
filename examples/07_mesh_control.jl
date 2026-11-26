# Mesh control: phase sizes, interface refinement with grading, a refined inclusion close to
# a periodic face (the size field is periodic), order-2 elements and the quality report.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))
rve = RVE(cell, [Phase(:stiff, Sphere((0.78, 0.5, 0.5), 0.2); refine = (h = 0.02, growth = 1.25)),
                 Phase(:soft, Sphere((0.25, 0.25, 0.25), 0.15); h = 0.04)])

res = generate(rve; h = 0.12, output = joinpath(@__DIR__, "mesh_control.msh"))
quality_report(res)

# quadratic elements with curved interfaces (for GridapGmsh with periodic cells use order 1)
res2 = generate(rve; mesh = MeshOptions(h = 0.12, order = 2))
quality_report(res2)
