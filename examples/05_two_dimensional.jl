# 2D cell: hexagonal packing of fibres seen in cross-section, plus an elliptical pore.
using RVEs

cell = Cell((1.0, sqrt(3.0)))
fibres = [Disk(cellpoint(cell, (0.5, 0.5)), 0.3), Disk(cellpoint(cell, (0.0, 0.0)), 0.3)]
pore = Rotate(Ellipse(cellpoint(cell, (0.5, 0.0)), (0.12, 0.05)); angle = π / 8)
rve = RVE(cell, [Phase(:fibre, fibres; refine = 0.02), Phase(:pore, pore)]; matrix = :resin)

res = generate(rve; h = 0.05, output = joinpath(@__DIR__, "hex2d.msh"))
display(res)
