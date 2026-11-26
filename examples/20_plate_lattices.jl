# Plate lattices (closed-cell): plates on the simple cubic, body-centred and face-centred cubic
# planes, and the elastically isotropic SC + FCC combination, by thickness or by relative density.
# The plates are what remains of the cell around the closed cells returned by plate_cells.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))
out(name) = joinpath(@__DIR__, name)
mesh = MeshOptions(h = 0.05, optimize = :netgen)

# single families with 20 % solid: the cells are empty (void) and the plates are the matrix
for kind in (:sc, :bcc, :fcc)
    t = plate_thickness(kind, 0.2)                       # thickness / edge
    cells = plate_cells(cell, kind; thickness = t)
    res = generate(RVE(cell, Phase(:cells, cells; void = true); matrix = :plates); mesh = mesh, output = out("plates_$kind.msh"))
    println(kind, " plates: thickness ", round(t; digits = 4), ", ", length(cells), " cells, density ",
            round(volume_fraction(res, :plates); digits = 4))
end

# isotropic SC + FCC ("cubic + octet", Berger et al. 2017): thicknesses SC : FCC = 8√3/9 : 1
t = plate_thickness((:sc, :fcc), 0.3)
res = generate(RVE(cell, Phase(:cells, plate_cells(cell, (:sc, :fcc); density = 0.3); void = true); matrix = :plates);
               mesh = mesh, output = out("plates_sc_fcc.msh"))
println("sc+fcc plates: thicknesses ", round.(t; digits = 4), ", density ", round(volume_fraction(res, :plates); digits = 4))

# SC + BCC with a chosen thickness ratio
t = plate_thickness((:sc, :bcc), 0.3; ratio = (1, 1))
println("sc+bcc plates with equal thicknesses for density 0.3: ", round.(t; digits = 4))

# own thicknesses, a non-cubic cell, and the cells filled with a second material
box = Cell((1.0, 1.0, 1.5); origin = (0.2, 0.0, -0.3))
cells = plate_cells(box, (:sc, :fcc); thickness = (0.08, 0.04))
res = generate(RVE(box, Phase(:filler, cells); matrix = :plates); mesh = mesh, output = out("plates_filled.msh"))
println("filled plate lattice: ", volume_fractions(res), "  exact plates ",
        round(plate_density(box, (:sc, :fcc); thickness = (0.08, 0.04)); digits = 4))
