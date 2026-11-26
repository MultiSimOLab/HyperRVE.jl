# Strut lattices: nodes in fractional coordinates, edges with integer cell offsets.

const _LTOL = 1e-9

_offsets3() = [(a, b, c) for a in -1:1 for b in -1:1 for c in -1:1]

"Minimum-image pairs of `nodes` at distance `d` (fractional, unit cube), each strut once."
function _edges_by_distance(nodes::Vector{Vec{3}}, d::Float64)
    edges = Tuple{Int,Int,NTuple{3,Int}}[]
    for i in eachindex(nodes), j in i:length(nodes), o in _offsets3()
        i == j && (o <= (0, 0, 0)) && continue          # each self-strut once, none of length 0
        dist = _norm(nodes[j] .+ o .- nodes[i])
        abs(dist - d) < 1e-9 && push!(edges, (i, j, o))
    end
    return edges
end

const _PRESET_NODES = Dict{Symbol,Vector{Vec{3}}}(
    :cubic => [(0.0, 0.0, 0.0)],
    :bcc => [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
    :fcc => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
    :octet => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
    :kelvin => unique([Tuple(p[k] for k in perm) for p in ((0.0, 0.25, 0.5), (0.0, 0.75, 0.5))
                       for perm in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))]),
    :diamond => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0),
                 (0.25, 0.25, 0.25), (0.25, 0.75, 0.75), (0.75, 0.25, 0.75), (0.75, 0.75, 0.25)])

const _PRESET_DIST = Dict(:cubic => 1.0, :bcc => sqrt(3) / 2, :fcc => sqrt(2) / 2, :octet => sqrt(2) / 2,
                          :kelvin => sqrt(2) / 4, :diamond => sqrt(3) / 4)

"""
    lattice_topology(kind) -> (nodes, edges)

Nodes (fractional coordinates) and struts `(i, j, offset)` of a preset strut lattice:
`:cubic`, `:bcc`, `:fcc` (corner-to-face-centre struts, i.e. the face diagonals, the usual
additive-manufacturing FCC), `:octet` (all nearest neighbours of FCC), `:kelvin`
(tetrakaidecahedral foam) and `:diamond`.
"""
function lattice_topology(kind::Symbol)
    haskey(_PRESET_NODES, kind) ||
        throw(ArgumentError("unknown lattice $(repr(kind)); presets: $(join(sort(collect(keys(_PRESET_NODES))), ", "))"))
    nodes = _PRESET_NODES[kind]
    edges = _edges_by_distance(nodes, _PRESET_DIST[kind])
    kind == :fcc && (edges = filter(e -> e[1] == 1 || e[2] == 1, edges))
    return (nodes=nodes, edges=edges)
end

_wrap01(x) = (y = x - floor(x); abs(y - 1) < 1e-12 ? 0.0 : y)

"""
    struts(cell, nodes, edges; radius, node_radius = :auto) -> Vector{Geometry{3}}

Items of a strut lattice: one [`Cylinder`](@ref) per strut and one [`Sphere`](@ref) per node.

- `nodes`: points in fractional coordinates of the cell (`(0.5, 0.5, 0.5)` is the centre).
- `edges`: `(i, j)` joins nodes `i` and `j` of the cell; `(i, j, (o1, o2, o3))` joins node `i`
  to the copy of node `j` in the cell shifted by `o` cells. List each strut once: periodic
  images are created by [`generate`](@ref).
- `radius`: one number or one per strut. `node_radius`: `:auto` (largest radius of the struts
  meeting at the node), a number, one per node, or `nothing` (no node spheres).

Along non-periodic directions every translate of a strut or node lying in the closed cell is
emitted, so the lattice is the complete finite truss. Use the items in a fused phase:
`Phase(:lattice, struts(...); fuse = true)`.
"""
function struts(cell::Cell{3}, nodes, edges; radius, node_radius=:auto)
    ξ = Vec{3}[_vec(Val(3), p, "lattice node") for p in nodes]
    isempty(ξ) && throw(ArgumentError("struts: no nodes"))
    # wrap nodes into [0, 1) and remember the shift to correct the edge offsets
    shift = [map(x -> floor(Int, x + 1e-12), p) for p in ξ]
    ξ = [map(_wrap01, p) for p in ξ]
    for i in eachindex(ξ), j in i+1:length(ξ), o in _offsets3()
        _norm(ξ[j] .+ o .- ξ[i]) < _LTOL && throw(ArgumentError("struts: nodes $i and $j coincide"))
    end
    E = Tuple{Int,Int,NTuple{3,Int}}[]
    for e in edges
        (length(e) in (2, 3)) || throw(ArgumentError("struts: an edge is (i, j) or (i, j, (o1, o2, o3)), got $e"))
        i, j = Int(e[1]), Int(e[2])
        o = length(e) == 3 ? ntuple(k -> Int(e[3][k]), 3) : (0, 0, 0)
        (1 <= i <= length(ξ) && 1 <= j <= length(ξ)) || throw(ArgumentError("struts: edge $e refers to a missing node"))
        o = o .+ Tuple(shift[j]) .- Tuple(shift[i])
        i == j && o == (0, 0, 0) && throw(ArgumentError("struts: edge $e has zero length"))
        push!(E, (i, j, o))
    end
    isempty(E) && throw(ArgumentError("struts: no edges"))
    canon(e) = (e[1], e[2], e[3]) <= (e[2], e[1], .-e[3]) ? (e[1], e[2], e[3]) : (e[2], e[1], .-e[3])
    seen = Set{Any}()
    for e in E
        c = canon(e)
        c in seen && throw(ArgumentError("struts: strut $(e) listed twice"))
        push!(seen, c)
    end
    r = radius isa Real ? fill(_positive(radius, "struts radius"), length(E)) :
        (length(radius) == length(E) ? [_positive(x, "struts radius") for x in radius] :
         throw(ArgumentError("struts: radius must be a number or one per strut ($(length(E)))")))
    2 * maximum(r) < minimum(cell.size) || throw(ArgumentError("struts: the radius must be smaller than half the cell"))
    # node radii
    rn = if node_radius === nothing
        nothing
    elseif node_radius === :auto
        v = zeros(length(ξ))
        for (k, (i, j, _)) in enumerate(E)
            v[i] = max(v[i], r[k]); v[j] = max(v[j], r[k])
        end
        v
    elseif node_radius isa Real
        fill(_positive(node_radius, "struts node_radius"), length(ξ))
    elseif length(node_radius) == length(ξ)
        [_positive(x, "struts node_radius") for x in node_radius]
    else
        throw(ArgumentError("struts: node_radius must be :auto, nothing, a number or one per node"))
    end
    used = falses(length(ξ))
    for (i, j, _) in E
        used[i] = used[j] = true
    end
    any(.!used) && @warn "struts: nodes $(findall(.!used)) have no struts and get no sphere"
    # translates along non-periodic directions (complete finite truss)
    per = cell.periodic
    inside(p) = all(d -> per[d] || (-_LTOL <= p[d] <= 1 + _LTOL), 1:3)
    svals = [per[d] ? (0,) : (-1, 0, 1) for d in 1:3]
    function placements(a, b)
        out = NTuple{3,Int}[]
        for s in Iterators.product(svals...)
            inside(a .+ s) && inside(b .+ s) && push!(out, s)
        end
        return isempty(out) ? [(0, 0, 0)] : out
    end
    items = Geometry{3}[]
    for (k, (i, j, o)) in enumerate(E)
        a, b = ξ[i], ξ[j] .+ o
        _norm(b .- a) <= 2 * r[k] && @warn "struts: strut $k is shorter than its diameter"
        for s in placements(a, b)
            p, q = cellpoint(cell, a .+ s), cellpoint(cell, b .+ s)
            push!(items, Cylinder(p, q .- p, r[k]))
        end
    end
    if rn !== nothing
        for i in eachindex(ξ)
            used[i] || continue
            for s in placements(ξ[i], ξ[i])
                push!(items, Sphere(cellpoint(cell, ξ[i] .+ s), rn[i]))
            end
        end
    end
    return items
end

# Relative density of the presets on a cubic cell, x = radius / edge: ρ = c1 x² − c2 x³,
# exact while the overlaps of different nodes do not touch (valid up to xmax); c2 fitted to
# OpenCASCADE volumes of the fused lattices (relative error below 5e-5).
const _DENSITY = Dict(:cubic => (3π, 8sqrt(2), 0.30), :bcc => (4sqrt(3) * π, 39.1918, 0.18),
                      :diamond => (4sqrt(3) * π, 39.192, 0.14), :kelvin => (6sqrt(2) * π, 60.334, 0.14),
                      :fcc => (6sqrt(2) * π, 54.6275, 0.14), :octet => (12sqrt(2) * π, 154.510, 0.14))

"""
    relative_density(kind, x)                         # preset, cubic cell, x = radius / edge
    relative_density(cell, nodes, edges; radius)      # any lattice, first order (upper bound)

Volume fraction of a strut lattice with node spheres of the strut radius. For the presets the
formula `c1 x² − c2 x³` is exact up to the radius where the overlaps of different nodes meet
(checked against OpenCASCADE volumes). For other lattices the first-order estimate
`Σ π r² L / V` overestimates by 10–20 % at practical radii: measure the exact value with
`volume_fraction(generate(rve; geometry_only = true), :lattice)`.
"""
function relative_density(kind::Symbol, x::Real)
    haskey(_DENSITY, kind) || throw(ArgumentError("unknown lattice $(repr(kind))"))
    c1, c2, xmax = _DENSITY[kind]
    0 < x <= xmax || throw(ArgumentError("relative_density(:$kind, x) is exact for 0 < x <= $xmax; " *
                                         "measure larger radii with generate(...; geometry_only = true)"))
    return c1 * x^2 - c2 * x^3
end
function relative_density(cell::Cell{3}, nodes, edges; radius)
    its = struts(cell, nodes, edges; radius=radius, node_radius=nothing)
    return sum(π * c.radius^2 * _norm(c.axis) for c in its if c isa Cylinder) / volume(cell)
end

"""
    strut_radius(kind, density) -> x

Strut radius divided by the cell edge that gives the relative `density` for a preset lattice
(inverse of [`relative_density`](@ref)).
"""
function strut_radius(kind::Symbol, ρ::Real)
    haskey(_DENSITY, kind) || throw(ArgumentError("unknown lattice $(repr(kind))"))
    c1, c2, xmax = _DENSITY[kind]
    0 < ρ <= relative_density(kind, xmax) ||
        throw(ArgumentError("density $ρ out of the range of :$kind (max $(round(relative_density(kind, xmax); digits=3)) " *
                            "with the exact formula); measure it with generate(...; geometry_only = true)"))
    x = sqrt(ρ / c1)
    for _ in 1:50
        f = c1 * x^2 - c2 * x^3 - ρ
        x -= f / (2c1 * x - 3c2 * x^2)
        abs(f) < 1e-15 && break
    end
    return x
end

"""
    lattice(cell, kind; radius = nothing, density = nothing, node_radius = :auto)

Items of a preset strut lattice (see [`lattice_topology`](@ref)) scaled to the cell, with
either the strut `radius` or the relative `density` (cubic cells only). Use them in a fused
phase: `Phase(:lattice, lattice(cell, :octet; density = 0.2); fuse = true)`.
"""
function lattice(cell::Cell{3}, kind::Symbol; radius=nothing, density=nothing, node_radius=:auto)
    (radius === nothing) == (density === nothing) && throw(ArgumentError("lattice: give exactly one of radius or density"))
    t = lattice_topology(kind)
    if density !== nothing
        maximum(cell.size) - minimum(cell.size) <= 1e-12 * maximum(cell.size) ||
            throw(ArgumentError("lattice: density needs a cubic cell; give the radius instead"))
        radius = strut_radius(kind, density) * cell.size[1]
    end
    return struts(cell, t.nodes, t.edges; radius=radius, node_radius=node_radius)
end
