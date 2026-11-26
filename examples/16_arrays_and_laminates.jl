# Benchmark microstructures: regular particle arrays and laminates.
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))
for kind in (:sc, :bcc, :fcc)
    res = generate(RVE(cell, Phase(:particles, particle_array(cell, kind; fraction = 0.3); refine = 0.03)); h = 0.1,
                   output = joinpath(@__DIR__, "array_$(kind).msh"))
    println(kind, ": ", length(res.phases[2].entities), " particle pieces, fraction ", round(volume_fraction(res, :particles); digits = 4))
end

# hexagonal array of fibres in 2D (cell 1 : √3)
c2 = Cell((1.0, sqrt(3)))
generate(RVE(c2, Phase(:fibres, particle_array(c2, :hex; fraction = 0.6)); matrix = :resin); h = 0.05,
         output = joinpath(@__DIR__, "array_hex2d.msh"))

# two-phase laminate normal to z, shifted so that no interface lies on a cell face
res = generate(RVE(cell, laminate(cell; axis = 3, fractions = (0.3, 0.7), names = (:stiff, :soft), offset = 0.1);
                   matrix = nothing); h = 0.1, output = joinpath(@__DIR__, "laminate.msh"))
println("laminate: ", [p.name => round(p.fraction; digits = 3) for p in res.phases])
