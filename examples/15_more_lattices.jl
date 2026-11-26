# BCCZ/FCCZ lattices and 2D strut lattices: square, triangular, honeycomb, kagome and the
# auxetic re-entrant honeycomb.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))
for kind in (:bccz, :fccz)
    res = generate(RVE(cell, Phase(:lattice, lattice(cell, kind; density = 0.12); fuse = true); matrix = nothing);
                   geometry_only = true)
    println(kind, ": strut radius ", round(strut_radius(kind, 0.12); digits = 4), ", density ",
            round(volume_fraction(res, :lattice); digits = 4))
end

# 2D lattices in their rectangular periodic cells (strut length a)
for kind in (:square, :triangular, :honeycomb, :kagome)
    c = lattice_cell(kind; a = 0.5)
    res = generate(RVE(c, Phase(:lattice, lattice(c, kind; thickness = 0.04); fuse = true)); h = 0.03,
                   output = joinpath(@__DIR__, "lattice2d_$(kind).msh"))
    println(kind, ": cell ", round.(c.size; digits = 3), ", solid fraction ", round(volume_fraction(res, :lattice); digits = 4))
end

# auxetic re-entrant honeycomb: inclined struts at 30° below the horizontal, vertical struts twice as long
c = lattice_cell(:reentrant; a = 0.5, angle = π / 6, ratio = 2.0)
res = generate(RVE(c, Phase(:lattice, lattice(c, :reentrant; thickness = 0.04, angle = π / 6, ratio = 2.0); fuse = true);
                   matrix = nothing); h = 0.02, output = joinpath(@__DIR__, "reentrant.msh"))
println("re-entrant: solid fraction ", round(volume_fraction(res, :lattice); digits = 4))
