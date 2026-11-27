# Triply periodic minimal surfaces (experimental): a gyroid sheet in a soft matrix and a
# Schwarz P network lattice. The TPMS is the only item of the only phase; its triangulated
# surface is meshed as it is.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))

# gyroid sheet at 30 % relative density in a matrix (interpenetrating composite)
g = tpms(:gyroid; type = :sheet, density = 0.30)
res = generate(RVE(cell, Phase(:gyroid, g)); output = joinpath(@__DIR__, "tpms_gyroid_sheet.msh"))
println(g, ": fraction ", round(volume_fraction(res, :gyroid); digits = 4),
        ", interface area ", round(interface_area(res, :matrix, :gyroid); digits = 3),
        ", min quality ", round(res.stats.min_quality; digits = 3))

# Schwarz P network lattice (no matrix), two periods per side
p = tpms(:schwarz_p; type = :network, density = 0.35, periods = 2)
res = generate(RVE(cell, Phase(:p, p); matrix = nothing); output = joinpath(@__DIR__, "tpms_p_network.msh"))
println(p, ": fraction ", round(volume_fraction(res, :p); digits = 4))

# levels for other densities
for ρ in (0.2, 0.3, 0.4)
    println("gyroid network level for density ", ρ, ": ", round(tpms_level(:gyroid; type = :network, density = ρ); digits = 4))
end
