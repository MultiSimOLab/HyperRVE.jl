# Auxetic 2D cells: chiral lattices (rings joined by tangent ligaments) and rotating rigid units.
using HyperRVE

# chiral and anti-chiral lattices: a = distance between the centres of connected rings
for kind in (:hexachiral, :tetrachiral, :trichiral, :antitetrachiral, :antitrichiral)
    c = lattice_cell(kind; a = 1.0)
    items = chiral(c, kind; ring_radius = 0.2, thickness = 0.04)
    res = generate(RVE(c, Phase(:lattice, items; fuse = true); matrix = nothing); h = 0.02,
                   output = joinpath(@__DIR__, "$(kind).msh"))
    println(kind, ": cell ", round.(c.size; digits = 3), ", solid fraction ", round(volume_fraction(res, :lattice); digits = 4),
            " (exact ", round(relative_density(kind; ring_radius = 0.2, thickness = 0.04); digits = 4), ")")
end

# the mirror image, with solid nodes, embedded in a soft matrix
c = lattice_cell(:hexachiral; a = 1.0)
items = chiral(c, :hexachiral; ring_radius = 0.2, thickness = 0.04, mirror = true, filled = true)
res = generate(RVE(c, Phase(:lattice, items; fuse = true)); h = 0.03, output = joinpath(@__DIR__, "hexachiral_in_matrix.msh"))
println("hexachiral in matrix: ", volume_fractions(res))

# rotating squares: side a, each square turned 30° from the closed state, hinges 0.08 wide
θ = π / 6
c = lattice_cell(:rotating_squares; a = 1.0, angle = θ)
res = generate(RVE(c, Phase(:squares, rotating_squares(c; angle = θ, hinge = 0.08); fuse = true); matrix = nothing);
               h = 0.04, output = joinpath(@__DIR__, "rotating_squares.msh"))
println("rotating squares: solid fraction ", round(volume_fraction(res, :squares); digits = 4),
        " (exact ", round(relative_density(:rotating_squares; angle = θ, hinge = 0.08); digits = 4), ")")

# rotating triangles
c = lattice_cell(:rotating_triangles; a = 1.0, angle = θ)
res = generate(RVE(c, Phase(:triangles, rotating_triangles(c; angle = θ, hinge = 0.08); fuse = true); matrix = nothing);
               h = 0.04, output = joinpath(@__DIR__, "rotating_triangles.msh"))
println("rotating triangles: solid fraction ", round(volume_fraction(res, :triangles); digits = 4),
        " (exact ", round(relative_density(:rotating_triangles; angle = θ, hinge = 0.08); digits = 4), ")")
