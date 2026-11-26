# Strut lattices: presets by relative density, a lattice embedded in a matrix, a custom
# topology and a finite (non-periodic) block.
using RVEs

cell = Cell((1.0, 1.0, 1.0))

# octet truss with 15 % solid, alone (porous) and inside a matrix
oct = Phase(:octet, lattice(cell, :octet; density = 0.15); fuse = true)
res = generate(RVE(cell, oct; matrix = nothing); h = 0.08, output = joinpath(@__DIR__, "octet.msh"))
println("octet: target 0.15, meshed ", round(volume_fraction(res, :octet); digits = 4))

res = generate(RVE(cell, Phase(:struts, lattice(cell, :bcc; radius = 0.07); fuse = true)); h = 0.08,
               output = joinpath(@__DIR__, "bcc_in_matrix.msh"))
println("bcc in matrix: ", interface_names(res))

# custom topology: nodes in fractional coordinates, struts (i, j, offset) to the image of node j
nodes = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]
edges = [(1, 2, (0, 0, 0)), (1, 2, (-1, 0, 0)), (1, 2, (0, -1, 0)), (1, 2, (-1, -1, 0)), (1, 1, (0, 0, 1))]
items = struts(cell, nodes, edges; radius = 0.06)
println("custom lattice, first-order density: ", round(relative_density(cell, nodes, edges; radius = 0.06); digits = 4))
res = generate(RVE(cell, Phase(:custom, items; fuse = true); matrix = nothing); h = 0.08)

# finite block of 2 x 1 x 1 BCC cells (no periodicity: the complete truss, with outer nodes)
blk = Cell((1.0, 1.0, 1.0); periodic = false)
res = generate(RVE(blk, Phase(:block, lattice(blk, :bcc; radius = 0.08); fuse = true); matrix = nothing,
                   tiling = (2, 1, 1)); h = 0.1, output = joinpath(@__DIR__, "bcc_block.msh"))
println("block volume: ", round(volume(res, :block); digits = 4))
