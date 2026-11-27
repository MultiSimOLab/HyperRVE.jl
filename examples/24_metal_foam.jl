# Closed-cell metal foam (Alporas-like, lengths in mm) with a pore size distribution:
# pore diameters from a truncated log-normal law, porosity 0.9, polyhedral and rounded pores.
using HyperRVE

dist = truncated_lognormal(4.0; cv = 0.35, min = 2.0, max = 8.0)     # pore diameter (mm)

# 30 pores in a cubic cell sized to hold them; polyhedral pores, wall from the porosity
f = laguerre_foam(; n = 30, diameter = dist, porosity = 0.9, rng = 1)
st = pore_statistics(f)
println(f, ", cell ", round(f.cell.size[1]; digits = 2), " mm")
println("pore diameter: mean ", round(st.mean; digits = 2), " (target ", round(st.target.mean; digits = 2), "), ",
        round(st.min; digits = 2), "-", round(st.max; digits = 2), " mm; ", round(st.faces_per_pore; digits = 1), " neighbours per pore")
res = generate(RVE(f.cell, Phase(:pores, pores(f); void = true); matrix = :walls); h = 1.0, retries = 3,
               output = joinpath(@__DIR__, "metal_foam.msh"))
println("polyhedral: solid fraction ", round(volume_fraction(res, :walls); digits = 4))

# a foam in the same cell with a 0.1 mm wall: the pores are rounded to reach the porosity
fr = laguerre_foam(f.cell; diameter = dist, porosity = 0.9, wall = 0.1, rng = 1)
res = generate(RVE(fr.cell, Phase(:pores, pores(fr); void = true); matrix = :walls); h = 1.0, retries = 3,
               output = joinpath(@__DIR__, "metal_foam_rounded.msh"))
println("rounded (", round(fr.rounding; digits = 3), "): solid fraction ", round(volume_fraction(res, :walls); digits = 4))
