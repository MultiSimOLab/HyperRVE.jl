# Unidirectional fibre composite: hexagonal packing of fibres along z, in a 3D cell.
using HyperRVE

a = 1.0
cell = Cell((a, sqrt(3) * a, 0.5a))
rf = 0.3a
fibres = [Fibre(cellpoint(cell, (0.5, 0.5, 0.0)), (0, 0, 1), rf),
          Fibre(cellpoint(cell, (0.0, 0.0, 0.0)), (0, 0, 1), rf)]   # corner fibre: 4 quarters
rve = RVE(cell, Phase(:fibre, fibres; refine = 0.04); matrix = :resin)

res = generate(rve; h = 0.12, output = joinpath(@__DIR__, "fibres.msh"))
display(res)
println("\nfibre volume fraction: ", volume_fraction(res, :fibre), " (exact ", 2π * rf^2 / (sqrt(3) * a^2), ")")
