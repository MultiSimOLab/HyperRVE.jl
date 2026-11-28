# Perforated metal foam: straight holes drilled through a closed-cell (Alporas-like) foam, the
# way such foams are made sound absorbing; one hole per cell or a lattice of them, with or
# without windows in the walls (lengths in mm).
using HyperRVE

dist = truncated_lognormal(4.0; cv = 0.35, min = 2.0, max = 8.0)
f = laguerre_foam(Cell((10.0, 10.0, 10.0)); diameter = dist, porosity = 0.9, wall = 0.1, rng = 1)

a = perforate(f; diameter = 1.5)                       # one hole per cell, centred, along z
r = a.report
println(r.holes, " hole, perforation ratio ", round(r.perforation_ratio; digits = 4), ", ", r.crossed_pores,
        " pores crossed, open porosity ", round(r.open_porosity; digits = 3), " of ", round(r.porosity; digits = 3))
res = generate(RVE(f.cell, foam_phases(a); matrix = :aluminium); h = 1.0, retries = 3,
               output = joinpath(@__DIR__, "perforated_foam.msh"))
println("air ", round(volume_fraction(res, :air); digits = 4), ", aluminium ", round(volume_fraction(res, :aluminium); digits = 4))

res = generate(RVE(f.cell, Phase(:air, a.fluid); matrix = nothing); h = 1.0, retries = 3,
               output = joinpath(@__DIR__, "perforated_foam_air.msh"))      # the fluid alone
println("air fraction ", round(volume_fraction(res, :air); digits = 4))

# 2 × 2 holes and windows in 20 % of the walls: one connected fluid network
b = perforate(f; diameter = 1.0, spacing = 5.0, windows = (open_fraction = 0.2,), rng = 3)
println(b.report.holes, " holes, ", b.report.windows, " windows, ", b.report.networks, " network, open porosity ",
        round(b.report.open_porosity; digits = 3))
res = generate(RVE(f.cell, foam_phases(b); matrix = :aluminium); h = 1.0, retries = 3,
               output = joinpath(@__DIR__, "perforated_foam_windows.msh"))
println("air ", round(volume_fraction(res, :air); digits = 4), ", aluminium ", round(volume_fraction(res, :aluminium); digits = 4))
