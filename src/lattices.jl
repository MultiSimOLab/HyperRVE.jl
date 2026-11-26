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
    :bccz => [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
    :fcc => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
    :fccz => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
    :octet => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
    :kelvin => unique([Tuple(p[k] for k in perm) for p in ((0.0, 0.25, 0.5), (0.0, 0.75, 0.5))
                       for perm in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))]),
    :diamond => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0),
                 (0.25, 0.25, 0.25), (0.25, 0.75, 0.75), (0.75, 0.25, 0.75), (0.75, 0.75, 0.25)])

const _PRESET_DIST = Dict(:cubic => 1.0, :bcc => sqrt(3) / 2, :fcc => sqrt(2) / 2, :bccz => sqrt(3) / 2, :fccz => sqrt(2) / 2, :octet => sqrt(2) / 2,
                          :kelvin => sqrt(2) / 4, :diamond => sqrt(3) / 4)

"""
    lattice_topology(kind) -> (nodes, edges)

Nodes (fractional coordinates) and struts `(i, j, offset)` of a preset strut lattice:
`:cubic`, `:bcc`, `:fcc` (corner-to-face-centre struts, i.e. the face diagonals, the usual
additive-manufacturing FCC), `:bccz` and `:fccz` (the same with struts along the z edges,
common in additive manufacturing), `:octet` (all nearest neighbours of FCC), `:kelvin`
(tetrakaidecahedral foam) and `:diamond`.
"""
function lattice_topology(kind::Symbol; kwargs...)
    kind in _PRESETS_2D && return lattice_topology(kind, Val(2); kwargs...)
    isempty(kwargs) || throw(ArgumentError("lattice_topology: angle and ratio apply to 2D lattices only"))
    haskey(_PRESET_NODES, kind) ||
        throw(ArgumentError("unknown lattice $(repr(kind)); presets: $(join(sort(collect(keys(_PRESET_NODES))), ", "))"))
    nodes = _PRESET_NODES[kind]
    edges = _edges_by_distance(nodes, _PRESET_DIST[kind])
    kind in (:fcc, :fccz) && (edges = filter(e -> e[1] == 1 || e[2] == 1, edges))
    kind in (:bccz, :fccz) && push!(edges, (1, 1, (0, 0, 1)))        # the z edges of the cube
    return (nodes=nodes, edges=edges)
end

_wrap01(x) = (y = x - floor(x); abs(y - 1) < 1e-12 ? 0.0 : y)

"""
    struts(cell, nodes, edges; radius, node_radius = :auto) -> Vector{Geometry}

Items of a strut lattice: one [`Cylinder`](@ref) per strut and one [`Sphere`](@ref) per node
(in a 2D cell: a rectangle of width `2 radius` per strut and a [`Disk`](@ref) per node;
`thickness = 2 radius` may be given instead of `radius`).

- `nodes`: points in fractional coordinates of the cell (`(0.5, 0.5, 0.5)` is the centre).
- `edges`: `(i, j)` joins nodes `i` and `j` of the cell; `(i, j, (o1, o2, o3))` joins node `i`
  to the copy of node `j` in the cell shifted by `o` cells. List each strut once: periodic
  images are created by [`generate`](@ref).
- `radius`: one number or one per strut. `node_radius`: `:auto` (largest radius of the struts
  meeting at the node), a number, one per node, or `nothing` (no node spheres).

Along non-periodic directions every translate of a strut or node that overlaps the closed
cell is emitted (a strut crossing a face appears at both faces, cut by them), so the lattice
is the complete finite truss. Use the items in a fused phase:
`Phase(:lattice, struts(...); fuse = true)`.
"""
function struts(cell::Cell{D}, nodes, edges; radius=nothing, thickness=nothing, node_radius=:auto) where {D}
    (radius === nothing) != (thickness === nothing) || throw(ArgumentError("struts: give radius or thickness (= 2 radius)"))
    (thickness === nothing || thickness isa Real || (thickness isa AbstractVector && all(x -> x isa Real, thickness))) ||
        throw(ArgumentError("struts: thickness must be a number or one per strut"))
    radius === nothing && (radius = thickness isa Real ? thickness / 2 : thickness ./ 2)
    ξ = Vec{D}[_vec(Val(D), p, "lattice node") for p in nodes]
    isempty(ξ) && throw(ArgumentError("struts: no nodes"))
    # wrap nodes into [0, 1) and remember the shift to correct the edge offsets
    shift = [map(x -> floor(Int, x + 1e-12), p) for p in ξ]
    ξ = [map(_wrap01, p) for p in ξ]
    for i in eachindex(ξ), j in i+1:length(ξ), o in Iterators.product(ntuple(_ -> -1:1, D)...)
        _norm(ξ[j] .+ o .- ξ[i]) < _LTOL && throw(ArgumentError("struts: nodes $i and $j coincide"))
    end
    E = Tuple{Int,Int,NTuple{D,Int}}[]
    for e in edges
        bad() = ArgumentError("struts: an edge is (i, j) or (i, j, offset) with $D integer offsets, got $e")
        (length(e) in (2, 3)) || throw(bad())
        all(x -> x isa Integer, (e[1], e[2])) || throw(bad())
        length(e) == 3 && !(length(e[3]) == D && all(x -> x isa Integer, e[3])) && throw(bad())
        i, j = Int(e[1]), Int(e[2])
        o = length(e) == 3 ? ntuple(k -> Int(e[3][k]), D) : ntuple(_ -> 0, D)
        (1 <= i <= length(ξ) && 1 <= j <= length(ξ)) || throw(ArgumentError("struts: edge $e refers to a missing node"))
        o = o .+ Tuple(shift[j]) .- Tuple(shift[i])
        i == j && all(==(0), o) && throw(ArgumentError("struts: edge $e has zero length"))
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
    radius isa Union{Real,AbstractVector} || throw(ArgumentError("struts: radius must be a number or one per strut"))
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
    elseif node_radius isa AbstractVector && length(node_radius) == length(ξ)
        [_positive(x, "struts node_radius") for x in node_radius]
    else
        throw(ArgumentError("struts: node_radius must be :auto, nothing, a number or one per node"))
    end
    used = falses(length(ξ))
    for (i, j, _) in E
        used[i] = used[j] = true
    end
    any(.!used) && @warn "struts: nodes $(findall(.!used)) have no struts and get no sphere"
    # translates along non-periodic directions (complete finite truss): along each such
    # direction, every translate whose extent overlaps the closed cell (or, for a strut parallel
    # to the faces, lies in it)
    per = cell.periodic
    function shifts(a, b, d)
        per[d] && return (0,)
        out = Int[]
        for s in -1:1
            lo, hi = minmax(a[d] + s, b[d] + s)
            ok = hi - lo < _LTOL ? (-_LTOL <= lo <= 1 + _LTOL) : (min(hi, 1.0) - max(lo, 0.0) > _LTOL)
            ok && push!(out, s)
        end
        return isempty(out) ? (0,) : Tuple(out)
    end
    placements(a, b) = vec(collect(Iterators.product(ntuple(d -> shifts(a, b, d), D)...)))
    items = Geometry{D}[]
    short = Int[]
    for (k, (i, j, o)) in enumerate(E)
        a, b = ξ[i], ξ[j] .+ o
        axis = (b .- a) .* cell.size                   # not q - p: exact lengths for the periodicity checks
        _norm(axis) <= 2 * r[k] && push!(short, k)
        for s in placements(a, b)
            p = cellpoint(cell, a .+ s)
            if D == 3
                push!(items, Cylinder(p, axis, r[k]))
            else                                       # a rectangle of width 2r along the strut
                n = (-axis[2], axis[1]) ./ _norm(axis) .* r[k]
                push!(items, Polygon([p .+ n, p .- n, p .+ axis .- n, p .+ axis .+ n]))
            end
        end
    end
    isempty(short) || @warn "struts: struts $short are shorter than their diameter"
    if rn !== nothing
        for i in eachindex(ξ)
            used[i] || continue
            for s in placements(ξ[i], ξ[i])
                push!(items, D == 3 ? Sphere(cellpoint(cell, ξ[i] .+ s), rn[i]) : Disk(cellpoint(cell, ξ[i] .+ s), rn[i]))
            end
        end
    end
    return items
end

# Relative density of the presets on a cubic cell, x = radius / edge: ρ = c1 x² − c2 x³,
# exact while the overlaps of different nodes do not touch (valid up to xmax); c2 fitted to
# OpenCASCADE volumes of the fused lattices (relative error below 5e-5).
const _DENSITY = Dict(:cubic => (3π, 8sqrt(2), 0.30), :bcc => (4sqrt(3) * π, 39.1918, 0.18),
                      :diamond => (4sqrt(3) * π, 39.192, 0.14), :kelvin => (6sqrt(2) * π, 60.338, 0.14),
                      :fcc => (6sqrt(2) * π, 54.6275, 0.14), :octet => (12sqrt(2) * π, 154.510, 0.14),
                      :bccz => ((4sqrt(3) + 1) * π, 49.1829, 0.18), :fccz => ((6sqrt(2) + 1) * π, 66.8563, 0.18))

"""
    relative_density(kind, x)                         # preset, cubic cell, x = radius / edge
                                                      # 2D presets: x = thickness / strut length
    relative_density(cell, nodes, edges; radius)      # any lattice, first order (upper bound)

Volume fraction of a strut lattice with node spheres of the strut radius. For the presets the
formula `c1 x² − c2 x³` is exact up to the radius where the overlaps of different nodes meet
(checked against OpenCASCADE volumes). For other lattices the first-order estimate
`Σ π r² L / V` overestimates by 10–20 % at practical radii (in non-periodic cells it counts
whole struts, also those the faces cut, and overestimates much more): measure the exact value with
`volume_fraction(generate(rve; geometry_only = true), :lattice)`.
"""
function relative_density(kind::Symbol, x::Real)
    if haskey(_DENSITY_2D, kind)
        0 < x < 0.5 || throw(ArgumentError("relative_density(:$kind, x): x = thickness / strut length must be in (0, 0.5)"))
        return _DENSITY_2D[kind](x)
    end
    haskey(_DENSITY, kind) || throw(ArgumentError("unknown lattice $(repr(kind))"))
    c1, c2, xmax = _DENSITY[kind]
    0 < x <= xmax * (1 + 1e-12) || throw(ArgumentError("relative_density(:$kind, x) is exact for 0 < x <= $xmax; " *
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
    return min(x, xmax)
end

"""
    lattice(cell, kind; radius = nothing, density = nothing, node_radius = :auto)

Items of a preset strut lattice (see [`lattice_topology`](@ref)) scaled to the cell, with
either the strut `radius` (a number or one per strut, in the order of `lattice_topology(kind).edges`)
or the relative `density` (cubic cells only; exact with the default `node_radius`, ranges up to
`relative_density(kind, xmax)`, see [`relative_density`](@ref)). Use them in a fused
phase: `Phase(:lattice, lattice(cell, :octet; density = 0.2); fuse = true)`.
"""
function lattice(cell::Cell{3}, kind::Symbol; radius=nothing, density=nothing, node_radius=:auto, kwargs...)
    isempty(kwargs) || throw(ArgumentError("lattice: unknown keywords $(join(keys(kwargs), ", ")) for a 3D lattice"))
    kind in _PRESETS_2D && throw(ArgumentError("lattice: $(repr(kind)) is a 2D lattice; use a 2D cell (lattice_cell($(repr(kind))))"))
    (radius === nothing) == (density === nothing) && throw(ArgumentError("lattice: give exactly one of radius or density"))
    t = lattice_topology(kind)
    if density !== nothing
        maximum(cell.size) - minimum(cell.size) <= 1e-12 * maximum(cell.size) ||
            throw(ArgumentError("lattice: density needs a cubic cell; give the radius instead"))
        radius = strut_radius(kind, density) * cell.size[1]
    end
    return struts(cell, t.nodes, t.edges; radius=radius, node_radius=node_radius)
end

# ---------------------------------------------------------------- 2D strut lattices

const _PRESETS_2D = (:square, :triangular, :honeycomb, :kagome, :reentrant)

"""
    lattice_topology(kind; angle, ratio) -> (nodes, edges, size)      # 2D presets (also lattice_topology(kind))

Nodes (fractional), struts `(i, j, offset)` and the rectangular periodic cell size, for a
strut length 1, of a 2D lattice: `:square`, `:triangular`, `:honeycomb` (hexagonal),
`:kagome` and `:reentrant` (auxetic re-entrant honeycomb: inclined struts of length 1 at
`angle` (radians) below the horizontal, default `π/6`, vertical struts of length `ratio`,
default 2, which must exceed `2 sin(angle)`;
`angle < 0` gives convex hexagons, `angle = -30°, ratio = 1` the regular honeycomb).
"""
function lattice_topology(kind::Symbol, ::Val{2}; angle=nothing, ratio=nothing)
    kind === :reentrant || (angle === nothing && ratio === nothing) ||
        throw(ArgumentError("lattice: angle and ratio apply to :reentrant only"))
    (angle === nothing || angle isa Real) && (ratio === nothing || ratio isa Real) ||
        throw(ArgumentError("lattice: angle (in radians) and ratio must be numbers"))
    angle = angle === nothing ? π / 6 : angle
    ratio = ratio === nothing ? 2.0 : ratio
    if kind === :square
        return (nodes=[(0.0, 0.0)], edges=[(1, 1, (1, 0)), (1, 1, (0, 1))], size=(1.0, 1.0))
    elseif kind === :triangular
        return (nodes=[(0.0, 0.0), (0.5, 0.5)],
                edges=[(1, 1, (1, 0)), (2, 2, (1, 0)), (1, 2, (0, 0)), (1, 2, (-1, 0)), (1, 2, (0, -1)), (1, 2, (-1, -1))],
                size=(1.0, sqrt(3)))
    elseif kind === :kagome
        L = (2.0, 2sqrt(3))
        P = [(0.0, 0.0), (1.0, 0.0), (0.5, sqrt(3) / 2), (1.0, sqrt(3)), (0.0, sqrt(3)), (1.5, 3sqrt(3) / 2)]
        ξ = [p ./ L for p in P]
        E = Tuple{Int,Int,NTuple{2,Int}}[]
        for i in eachindex(ξ), j in i:length(ξ), o in Iterators.product(-1:1, -1:1)
            i == j && o <= (0, 0) && continue
            abs(_norm((ξ[j] .+ o .- ξ[i]) .* L) - 1) < 1e-9 && push!(E, (i, j, o))
        end
        return (nodes=ξ, edges=E, size=L)
    elseif kind in (:honeycomb, :reentrant)
        θ, h = kind === :honeycomb ? (-π / 6, 1.0) : (Float64(angle), Float64(ratio))
        (h > 0 && abs(θ) < π / 2) || throw(ArgumentError("lattice: ratio must be > 0 and |angle| < 90°"))
        W, H = 2cos(θ), 2h - 2sin(θ)
        h > 2max(sin(θ), 0.0) || throw(ArgumentError("lattice: for a re-entrant angle the ratio must exceed 2 sin(angle) = " *
                                                    "$(round(2sin(θ); digits=4)) (the re-entrant tips would touch)"))
        nodes = [(0.0, 0.0), (0.0, h / H), (0.5, (h - sin(θ)) / H), (0.5, (2h - sin(θ)) / H)]
        edges = [(1, 2, (0, 0)), (3, 4, (0, 0)), (2, 3, (0, 0)), (2, 3, (-1, 0)), (4, 1, (0, 1)), (4, 1, (1, 1))]
        # nodes in [0, 1): wrap and correct the offsets of the edges
        sh = [floor.(Int, p .+ 1e-12) for p in nodes]
        nodes = [p .- q for (p, q) in zip(nodes, sh)]
        edges = [(i, j, Tuple(o .+ sh[j] .- sh[i])) for (i, j, o) in edges]
        return (nodes=nodes, edges=edges, size=(W, H))
    end
    throw(ArgumentError("unknown 2D lattice $(repr(kind)); presets: $(join(_PRESETS_2D, ", "))"))
end

"""
    lattice_cell(kind; a = 1, angle, ratio, origin = nothing) -> Cell{2}

The rectangular periodic cell of a 2D preset lattice with strut length `a` (see
[`lattice_topology`](@ref)); `lattice(lattice_cell(:honeycomb; a = 0.2), :honeycomb; thickness = 0.02)`.
"""
function lattice_cell(kind::Symbol; a=1.0, angle=nothing, ratio=nothing, origin=nothing)
    t = lattice_topology(kind, Val(2); angle, ratio)
    return Cell(_positive(a, "lattice_cell a") .* t.size; origin=origin)
end

"""
    lattice(cell::Cell{2}, kind; thickness, node_radius = :auto, angle, ratio)

Items of a 2D preset strut lattice (rectangles of width `thickness` and disks at the nodes)
in its rectangular periodic cell, which must have the proportions of
[`lattice_cell`](@ref)`(kind; ...)` (any size). Use them in a fused phase:
`Phase(:lattice, lattice(cell, :reentrant; thickness = 0.05); fuse = true)`.
"""
function lattice(cell::Cell{2}, kind::Symbol; thickness=nothing, radius=nothing, node_radius=:auto, angle=nothing, ratio=nothing)
    (thickness === nothing) != (radius === nothing) || throw(ArgumentError("lattice: give thickness (or radius = thickness / 2)"))
    t = lattice_topology(kind, Val(2); angle, ratio)
    abs(cell.size[1] / cell.size[2] - t.size[1] / t.size[2]) <= 1e-6 * t.size[1] / t.size[2] ||
        throw(ArgumentError("lattice: the cell of a $(repr(kind)) lattice must have the proportions " *
                            "$(round(t.size[1] / t.size[2]; sigdigits=6)) (width / height); use lattice_cell($(repr(kind)); a = ...)"))
    return struts(cell, t.nodes, t.edges; thickness=thickness, radius=radius, node_radius=node_radius)
end

# Area fraction of the 2D presets with strut length a and thickness t, x = t / a (exact, with or
# without the node disks, while the struts of a node do not touch other struts)
const _DENSITY_2D = Dict(:square => x -> 2x - x^2,
                         :triangular => x -> 1 - (1 - sqrt(3) * x)^2,
                         :honeycomb => x -> 1 - (1 - x / sqrt(3))^2,
                         :kagome => x -> 1 - 0.75 * (1 - x / sqrt(3))^2 - 0.25 * (1 - sqrt(3) * x)^2)
