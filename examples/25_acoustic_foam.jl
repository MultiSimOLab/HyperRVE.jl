# Acoustic metal foam: windows in a fraction of the walls between pores and the connected
# fluid (air) phase, for cell problems of porous acoustics (lengths in mm).
using HyperRVE

dist = truncated_lognormal(4.0; cv = 0.35, min = 2.0, max = 8.0)
f = laguerre_foam(; n = 16, diameter = dist, porosity = 0.9, wall = 0.1, rng = 1)   # rounded pores
a = open_pores(f; open_fraction = 0.3, window = 0.4, rng = 3)    # 30 % of the walls get a window
r = a.report
println(r.windows, " windows (", round(r.open_area_fraction; digits = 3), " of the wall area), ",
        round(r.windows_per_pore; digits = 1), " per pore; percolates ", r.percolates,
        "; open porosity ", round(r.open_porosity; digits = 4), ", ", r.closed_pores, " closed pores")

# the fluid alone (its wall is the group interface_air_void)
res = generate(RVE(f.cell, Phase(:air, a.fluid); matrix = nothing); h = 1.0, retries = 3,
               output = joinpath(@__DIR__, "acoustic_foam_air.msh"))
println("air fraction ", round(volume_fraction(res, :air); digits = 4), ", groups ", filter(n -> startswith(n, "interface"), interface_names(res)))

# fluid and solid (closed pores, if any, stay empty)
phases = [Phase(:air, a.fluid)]
isempty(a.closed) || push!(phases, Phase(:closed, a.closed; void = true))
res = generate(RVE(f.cell, phases; matrix = :aluminium); h = 1.0, retries = 3,
               output = joinpath(@__DIR__, "acoustic_foam.msh"))
println("air ", round(volume_fraction(res, :air); digits = 4), ", aluminium ", round(volume_fraction(res, :aluminium); digits = 4))
