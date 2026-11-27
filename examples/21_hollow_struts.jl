# Hollow strut lattices (tubes): the lattice of radius R painted over by its cavity, a lattice of
# radius R − wall (hollow nodes, one connected cavity) or one cylinder per strut (solid nodes).
using HyperRVE

cell = Cell((1.0, 1.0, 1.0))
out(name) = joinpath(@__DIR__, name)

# hollow nodes, empty inside: tubes of radius 0.1 and wall 0.03
hl = hollow_lattice(cell, :cubic; radius = 0.1, wall = 0.03)
rve = RVE(cell, [Phase(:lattice, hl.shell; fuse = true), Phase(:cavity, hl.core; fuse = true, void = true)]; matrix = nothing)
res = generate(rve; h = 0.04, output = out("hollow_cubic.msh"))
println("hollow cubic: density ", round(volume_fraction(res, :lattice); digits = 4), " (solid struts: ",
        round(relative_density(:cubic, 0.1); digits = 4), ")")

# solid nodes: the cavities are closed, one per strut; wall from the relative density
hs = hollow_lattice(cell, :bcc; radius = 0.08, density = 0.07, node = :solid)
rve = RVE(cell, [Phase(:lattice, hs.shell; fuse = true), Phase(:cavity, hs.core; void = true)]; matrix = nothing)
res = generate(rve; h = 0.04, output = out("hollow_bcc_solid_nodes.msh"))
println("bcc, solid nodes: target 0.07, meshed ", round(volume_fraction(res, :lattice); digits = 4))

# tubes in a matrix with the inside filled by a third material (a core-shell lattice)
rve = RVE(cell, [Phase(:shell, hl.shell; fuse = true), Phase(:core, hl.core; fuse = true)])
res = generate(rve; h = 0.05, output = out("core_shell_cubic.msh"))
println("core-shell: ", volume_fractions(res))
println(interface_names(res))
