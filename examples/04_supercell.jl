# A 3 x 3 x 3 supercell of a unit cell containing a rotated ellipsoid.
using RVEs

cell = Cell((1.0, 1.0, 1.0))
e = Rotate(Ellipsoid(cellpoint(cell, (0.5, 0.5, 0.5)), (0.35, 0.12, 0.12)); axis = (1, 1, 0), angle = π / 4)
unit = RVE(cell, Phase(:inclusion, e))

res = generate(supercell(unit, (3, 3, 3)); h = 0.15, output = joinpath(@__DIR__, "supercell.msh"))
display(res)
