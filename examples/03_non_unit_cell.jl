# Cells of any size and origin: a 2 x 1 x 0.5 mm cell placed at (10, -2, 3) mm, periodic in x
# and y only, with a rotated ellipsoid crossing a periodic face and a sphere cut by the free z face.
using RVEs

cell = Cell((2.0, 1.0, 0.5); origin = (10.0, -2.0, 3.0), periodic = (true, true, false))
e = Rotate(Ellipsoid(cellpoint(cell, (0.95, 0.5, 0.5)), (0.4, 0.15, 0.15)); axis = (0, 0, 1), angle = π / 6)   # crosses xmax
s = Sphere(cellpoint(cell, (0.5, 0.5, 1.0)), 0.2)             # cut by the non-periodic zmax face
rve = RVE(cell, [Phase(:ellipsoid, e; refine = 0.03), Phase(:sphere, s; refine = 0.03)])

res = generate(rve; h = 0.1, output = joinpath(@__DIR__, "non_unit.msh"))
display(res)
