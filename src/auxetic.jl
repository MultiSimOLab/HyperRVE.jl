# ---------------------------------------------------------------- auxetic 2D cells

const _CHIRAL = (:hexachiral, :tetrachiral, :trichiral, :antitetrachiral, :antitrichiral)
const _ROTATING = (:rotating_squares, :rotating_triangles)
const _AUXETIC = (_CHIRAL..., _ROTATING...)

# ring centres (fractions of the cell), connections, cell size for unit ring spacing and
# ring colours (nothing for chiral; ±1 for anti-chiral, the side of the ligaments)
function _chiral_topology(kind::Symbol)
    if kind === :hexachiral
        t = lattice_topology(:triangular, Val(2))
        return (t.nodes, t.edges, t.size, nothing)
    elseif kind === :tetrachiral
        t = lattice_topology(:square, Val(2))
        return (t.nodes, t.edges, t.size, nothing)
    elseif kind === :trichiral || kind === :antitrichiral
        t = lattice_topology(:honeycomb, Val(2))
        return (t.nodes, t.edges, t.size, kind === :trichiral ? nothing : [1, -1, 1, -1])
    elseif kind === :antitetrachiral
        nodes = [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)]
        edges = [(1, 2, (0, 0)), (2, 1, (1, 0)), (3, 4, (0, 0)), (4, 3, (1, 0)),
                 (1, 3, (0, 0)), (3, 1, (0, 1)), (2, 4, (0, 0)), (4, 2, (0, 1))]
        return (nodes, edges, (2.0, 2.0), [1, -1, -1, 1])
    end
    throw(ArgumentError("unknown chiral lattice $(repr(kind)); kinds: $(join(_CHIRAL, ", ")) " *
                        "(rotating units: $(join(_ROTATING, ", ")))"))
end

_rot_scale(kind, θ) = kind === :rotating_squares ? 2 * (cos(θ) + sin(θ)) : cos(θ) + sqrt(3) * sin(θ)

function _rot_angle(kind, angle)
    angle isa Real || throw(ArgumentError("$kind: give the rotation angle in radians"))
    θmax = kind === :rotating_squares ? π / 2 : 2π / 3
    0 < angle < θmax || throw(ArgumentError("$kind: angle must be in (0, $(kind === :rotating_squares ? "π/2" : "2π/3"))"))
    return float(angle)
end

# cell size for unit ring spacing / unit side, used by lattice_cell
function _auxetic_size(kind, angle)
    kind in _CHIRAL && return _chiral_topology(kind)[3]
    θ = _rot_angle(kind, angle)
    s = _rot_scale(kind, θ)
    return kind === :rotating_squares ? (s, s) : (s, sqrt(3) * s)
end

function _auxetic_cell(cell::Cell{2}, kind, sz)
    all(cell.periodic) || throw(ArgumentError("$kind: the cell must be periodic in both directions; " *
                                              "build finite samples with supercell"))
    abs(cell.size[1] / cell.size[2] - sz[1] / sz[2]) <= 1e-6 * sz[1] / sz[2] ||
        throw(ArgumentError("$kind: the cell must have the proportions $(round(sz[1] / sz[2]; sigdigits=6)) " *
                            "(width / height); use lattice_cell($(repr(kind)); a = ...$(kind in _ROTATING ? ", angle = ..." : ""))"))
    return cell.size[1] / sz[1]
end

function _check_chiral(kind, r, t, filled, a)
    (r isa Real && t isa Real) || throw(ArgumentError("$kind: give ring_radius and thickness"))
    r > 0 && t > 0 || throw(ArgumentError("$kind: ring_radius and thickness must be positive"))
    filled || t < 2r || throw(ArgumentError("$kind: thickness must be smaller than 2 ring_radius (or use filled = true)"))
    r + t / 2 < a / 2 || throw(ArgumentError("$kind: neighbouring rings touch; need ring_radius + thickness / 2 < " *
                                            "$(a / 2) (half the ring spacing)"))
    return nothing
end

"""
    chiral(cell::Cell{2}, kind; ring_radius, thickness, mirror = false, filled = false)

Items of a chiral or anti-chiral lattice: rings of mean radius `ring_radius` and wall
`thickness`, joined by straight ligaments of the same thickness tangent to the rings
(outer radius `ring_radius + thickness / 2`). `kind` is `:hexachiral`, `:tetrachiral`,
`:trichiral`, `:antitetrachiral` or `:antitrichiral`; the cell must be periodic with the
proportions of [`lattice_cell`](@ref)`(kind)`, whose `a` is the distance between the
centres of connected rings. `mirror = true` gives the mirror image (opposite handedness);
the anti-chiral lattices are their own mirror image, so for them it only moves the lattice.
`filled = true` gives solid disks instead of rings. Use them in a fused phase:
`Phase(:lattice, chiral(c, :hexachiral; ring_radius = 0.2, thickness = 0.04); fuse = true)`.
The exact solid fraction is [`relative_density`](@ref)`(kind; ring_radius, thickness, a)`.
"""
function chiral(cell::Cell{2}, kind::Symbol; ring_radius=nothing, thickness=nothing, mirror::Bool=false, filled::Bool=false)
    nodes, edges, sz, col = _chiral_topology(kind)
    a = _auxetic_cell(cell, kind, sz)
    r, t = ring_radius, thickness
    _check_chiral(kind, r, t, filled, a)
    # anti-chiral ligaments run along the rows of rings, which lie on cell faces: with t = 2r a
    # ligament edge would lie on the face and break the periodic matching
    col === nothing || abs(t - 2r) > 1e-6 * a ||
        throw(ArgumentError("$kind: thickness = 2 ring_radius puts the ligament edges on the cell faces; " *
                            "change the thickness slightly"))
    s0 = mirror ? -1 : 1
    items = Geometry{2}[]
    for ξ in nodes
        c = cellpoint(cell, ξ)
        push!(items, filled ? Disk(c, r + t / 2) : Cut(Disk(c, r + t / 2), Disk(c, r - t / 2)))
    end
    for (i, j, o) in edges
        A = cellpoint(cell, nodes[i])
        v = (nodes[j] .+ o .- nodes[i]) .* cell.size
        R = _norm(v)
        d = v ./ R
        n = (-d[2], d[1])
        sβ = col === nothing ? 2r / R : 0.0            # chiral: ligament leaves the ring at the tangent point
        s = col === nothing ? s0 : s0 * col[i]
        m = sβ .* d .+ (s * sqrt(1 - sβ^2)) .* n
        P = A .+ r .* m
        Q = col === nothing ? A .+ v .- r .* m : A .+ v .+ r .* m
        u = (Q .- P) ./ _norm(Q .- P)
        w = (-u[2], u[1]) .* (t / 2)
        push!(items, Polygon([P .+ w, P .- w, Q .- w, Q .+ w]))
    end
    return items
end

_ngon(c, a, n, φ0) = [c .+ (a / (2sin(π / n))) .* (cos(φ0 + 2π * k / n), sin(φ0 + 2π * k / n)) for k in 0:n-1]

function _check_hinge(kind, hinge, θ, a)
    hinge isa Real || throw(ArgumentError("$kind: give the hinge width"))
    lim = kind === :rotating_squares ? 2a * min(sin(θ), cos(θ)) : 2a * min(sin(θ), sin(2π / 3 - θ))
    0 < hinge < min(a, lim) || throw(ArgumentError("$kind: hinge must be in (0, $(min(a, lim))) for this angle " *
                                                   "(wider hinges merge neighbouring units)"))
    return nothing
end

"""
    rotating_squares(cell::Cell{2}; angle, hinge)

Items of the rotating-squares mechanism: squares turned by `±angle` from the closed state
(alternating), joined at their corners by hinge disks of diameter `hinge`. The cell must be
periodic and square, [`lattice_cell`](@ref)`(:rotating_squares; a, angle)`, with `a` the side
of the squares; `a` follows from the cell size and `angle`, so give the same `angle` to both.
`angle` is in radians, in (0, π/2); the openings between squares have angle `2 angle`, and `angle`
and `π/2 - angle` give mirror-image cells. Exact solid fraction:
[`relative_density`](@ref)`(:rotating_squares; angle, hinge, a)`.
"""
function rotating_squares(cell::Cell{2}; angle=nothing, hinge=nothing)
    θ = _rot_angle(:rotating_squares, angle)
    s = _rot_scale(:rotating_squares, θ)
    a = _auxetic_cell(cell, :rotating_squares, (s, s))
    _check_hinge(:rotating_squares, hinge, θ, a)
    d = cell.size[1] / 2
    items = Geometry{2}[]
    for (i, j) in ((0, 0), (1, 0), (0, 1), (1, 1))
        sg = iseven(i + j) ? 1 : -1
        v = _ngon(cell.origin .+ (i * d, j * d), a, 4, π / 4 + sg * θ)
        push!(items, Polygon(v))
        sg == 1 && append!(items, [Disk(p, hinge / 2) for p in v])
    end
    return items
end

"""
    rotating_triangles(cell::Cell{2}; angle, hinge)

Items of the rotating-triangles mechanism: equilateral triangles of side `a` turned by
`angle` from the closed (triangular tiling) state, joined at their corners by hinge disks of
diameter `hinge`. The cell must be periodic with the proportions of
[`lattice_cell`](@ref)`(:rotating_triangles; a, angle)`; `a` follows from the cell size and
`angle`, so give the same `angle` to both. `angle` is in radians, in (0, 2π/3); `π/3` is the
kagome arrangement. Exact solid fraction: [`relative_density`](@ref)`(:rotating_triangles; angle, hinge, a)`.
"""
function rotating_triangles(cell::Cell{2}; angle=nothing, hinge=nothing)
    θ = _rot_angle(:rotating_triangles, angle)
    s = _rot_scale(:rotating_triangles, θ)
    a = _auxetic_cell(cell, :rotating_triangles, (s, sqrt(3) * s))
    _check_hinge(:rotating_triangles, hinge, θ, a)
    b = cell.size[1]
    φ = π / 3 - θ
    items = Geometry{2}[]
    for P in ((0.0, 0.0), (b / 2, sqrt(3) * b / 2))
        v = _ngon(cell.origin .+ P, a, 3, π / 2 + φ)
        push!(items, Polygon(v))
        append!(items, [Disk(p, hinge / 2) for p in v])
        push!(items, Polygon(_ngon(cell.origin .+ P .+ (b / 2, -b / (2sqrt(3))), a, 3, -π / 2 - φ)))
    end
    return items
end

"""
    relative_density(kind; ring_radius, thickness, filled = false, a = 1)   # chiral kinds
    relative_density(kind; angle, hinge, a = 1)                             # :rotating_squares, :rotating_triangles

Exact solid area fraction of the auxetic cells, with `a` the ring spacing or the side of the
units (the same lengths as in [`chiral`](@ref), [`rotating_squares`](@ref) and
[`rotating_triangles`](@ref)). It is exact while neighbouring ligaments do not overlap:
`thickness < 2 ring_radius / 3` for `:hexachiral`, and with `filled = true` also
`thickness <= 2 ring_radius` (tetra-) or `thickness <= 6 ring_radius` (tri-); otherwise it throws.
Measure other cases with `volume_fraction(generate(rve; geometry_only = true), name)`.
"""
function relative_density(kind::Symbol; ring_radius=nothing, thickness=nothing, filled::Bool=false, angle=nothing, hinge=nothing, a::Real=1.0)
    a > 0 || throw(ArgumentError("relative_density: a must be positive"))
    if kind in _ROTATING
        (ring_radius === nothing && thickness === nothing && !filled) ||
            throw(ArgumentError("relative_density(:$kind): give angle and hinge (ring_radius, thickness, filled are for chiral lattices)"))
        θ = _rot_angle(kind, angle)
        _check_hinge(kind, hinge, θ, a)
        w = hinge / a
        s = _rot_scale(kind, θ)
        kind === :rotating_squares && return (1 + π * w^2 / 4) / (s / 2)^2
        return (sqrt(3) + π * w^2) / (sqrt(3) * s^2)
    end
    kind in _CHIRAL || throw(ArgumentError("relative_density: unknown auxetic cell $(repr(kind)); kinds: $(join(_AUXETIC, ", ")) " *
                                           "(strut lattices: relative_density(kind, x))"))
    (angle === nothing && hinge === nothing) ||
        throw(ArgumentError("relative_density(:$kind): give ring_radius and thickness (angle, hinge are for rotating units)"))
    nodes, edges, sz, col = _chiral_topology(kind)
    _check_chiral(kind, ring_radius, thickness, filled, a)
    r, t = ring_radius / a, thickness / a
    limit = kind === :hexachiral ? (2r / 3, "<", "2 ring_radius / 3") :
            !filled ? (Inf, "<", "") :
            kind in (:tetrachiral, :antitetrachiral) ? (2r, "<=", "2 ring_radius") : (6r, "<=", "6 ring_radius")
    (limit[2] == "<" ? t < limit[1] : t <= limit[1] * (1 + 1e-12)) ||
        throw(ArgumentError("relative_density(:$kind$(filled ? "; filled = true" : "")): exact for thickness $(limit[2]) $(limit[3]), " *
                            "where neighbouring ligaments do not overlap; measure thicker walls with generate(...; geometry_only = true)"))
    L = col === nothing ? sqrt(1 - 4r^2) : 1.0
    Ro, Ri = r + t / 2, r - t / 2
    Aov = π / 4 * Ro^2 - Ri / 2 * sqrt(Ro^2 - Ri^2) - Ro^2 / 2 * asin(Ri / Ro)
    ring = filled ? π * Ro^2 : 2π * r * t
    return (length(nodes) * ring + length(edges) * (L * t - 2Aov)) / prod(sz)
end
