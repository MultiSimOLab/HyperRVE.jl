# Foams from a Voronoi tessellation: open-cell (edges as struts) and closed-cell (shrunk grains
# as pores, walls of constant thickness).

"""
    tessellation_edges(t::Tessellation{3}) -> (nodes, edges)

Vertices (fractional) and edges `(i, j, offset)` of the tessellation, each edge once (the
input of [`struts`](@ref)).
"""
function tessellation_edges(t::Tessellation{3})
    E = Set{Tuple{Int,Int,NTuple{3,Int}}}()
    for loops in t.loops, f in loops, i in eachindex(f)
        (v1, o1), (v2, o2) = f[i], f[mod1(i + 1, length(f))]
        e, r = (v1, v2, o2 .- o1), (v2, v1, o1 .- o2)
        push!(E, min(e, r))
    end
    return t.vertices, sort!(collect(E))
end

"""
    foam_struts(t; radius, node_radius = :auto, min_length = nothing) -> Vector{Geometry{3}}

Open-cell foam: the edges of a fully periodic tessellation as struts (one `Cylinder` per edge
and one `Sphere` per vertex, see [`struts`](@ref)). Edges shorter than `min_length` are
collapsed first (their vertices merged; `min_length = nothing` means `2 radius`), since a strut shorter than its diameter is buried in
the node spheres. Use in a fused phase without matrix:
`RVE(cell, Phase(:foam, foam_struts(t; radius = 0.02); fuse = true); matrix = nothing)`.
Practical for about 5–20 cells (the fused union of the struts dominates the cost).
"""
function foam_struts(t::Tessellation{3}; radius, node_radius=:auto, min_length=nothing)
    all(t.periodic) || throw(ArgumentError("foam_struts: the tessellation must be periodic in every direction " *
                                           "(build it with voronoi(cell; ..., periodic = true))"))
    r = _positive(radius, "foam_struts radius")
    lmin = min_length === nothing ? 2r : Float64(min_length)
    lmin > 0.25 * minimum(t.cell.size) &&
        @warn "foam_struts: min_length = $lmin is large compared with the cell; most vertices will be merged"
    nodes, edges = tessellation_edges(t)
    L = t.cell.size
    # union–find over vertices; T(v) is the lattice offset of the copy of v that belongs to the
    # cluster at its root's position (tp: the same, relative to the parent)
    parent = collect(eachindex(nodes))
    tp = [(0, 0, 0) for _ in nodes]
    function find(v)
        parent[v] == v && return v, (0, 0, 0)
        root, Tp = find(parent[v])
        tp[v] = tp[v] .+ Tp
        parent[v] = root
        return root, tp[v]
    end
    for (a, b, o) in edges
        _norm((nodes[b] .+ o .- nodes[a]) .* L) < lmin || continue
        (ra, Ta), (rb, Tb) = find(a), find(b)
        ra == rb && continue
        parent[rb] = ra
        tp[rb] = o .+ Ta .- Tb
    end
    T = [find(v)[2] for v in eachindex(nodes)]
    roots = sort!(unique(parent))
    id = Dict(q => i for (i, q) in enumerate(roots))
    pos = [(0.0, 0.0, 0.0) for _ in roots]
    cnt = zeros(Int, length(roots))
    for v in eachindex(nodes)
        i = id[parent[v]]
        pos[i] = pos[i] .+ nodes[v] .+ T[v]; cnt[i] += 1
    end
    pos = [pos[i] ./ cnt[i] for i in eachindex(pos)]
    out = Set{Tuple{Int,Int,NTuple{3,Int}}}()
    for (a, b, o) in edges
        ia, ib = id[parent[a]], id[parent[b]]
        k = o .- T[b] .+ T[a]
        ia == ib && all(==(0), k) && continue           # collapsed
        push!(out, min((ia, ib, k), (ib, ia, .-k)))
    end
    return struts(t.cell, pos, sort!(collect(out)); radius=r, node_radius=node_radius)
end

# ---------------------------------------------------------------- closed-cell foams

# grain i as a clipping polyhedron/polygon in cell-local coordinates, with the outward plane of
# each face (normal, offset)
function _grain_planes(t::Tessellation{3}, i)
    v, f = _grain_geometry(t, i)
    s = _polyhedron_centroid(v, f)                     # inside (a Laguerre seed may lie outside its grain)
    planes = Tuple{NTuple{3,Float64},Float64}[]
    for face in f
        n = _face_normal(v, face)
        n = n ./ _norm(n)
        c = sum(n .* v[face[1]])
        sum(n .* s) > c && (n = .-n; c = -c)            # outward: the centroid is on the inner side
        push!(planes, (n, c))
    end
    return _Poly3([[v[k] for k in face] for face in f], [(0, (0, 0, 0)) for _ in f]), planes
end
function _grain_planes(t::Tessellation{2}, i)
    loop = only(t.loops[i])
    v = [_vpos(t, vo) for vo in loop]
    s = _polygon_centroid(v)
    planes = Tuple{NTuple{2,Float64},Float64}[]
    for k in eachindex(v)
        p, q = v[k], v[mod1(k + 1, length(v))]
        n = (q[2] - p[2], p[1] - q[1]); n = n ./ _norm(n)
        c = sum(n .* p)
        sum(n .* s) > c && (n = .-n; c = -c)
        push!(planes, (n, c))
    end
    return _Poly2(v, [(0, (0, 0)) for _ in v]), planes
end

function _shrunk(t::Tessellation{D}, i, w) where {D}
    P, planes = _grain_planes(t, i)
    w == 0 && return P
    eps = 1e-11 * maximum(t.cell.size)
    for (n, c) in planes
        _clip!(P, n, c - w, (0, ntuple(_ -> 0, D)), eps)
        (D == 3 ? isempty(P.faces) : length(P.pts) < 3) && return nothing
    end
    return P
end

# the weld tolerance is relative to the size of the body (cell-local coordinates, about the cell size)
function _poly_volume(P::_Poly3)
    sc = maximum(p -> maximum(abs, p), _vertices(P))
    return abs(_polyhedron_volume(_indexed(P, 1e-12 * sc)...))
end
_poly_volume(P::_Poly2) = abs(_polygon_area(P.pts))

# pieces of a body (cell-local clipping polyhedron/polygon) inside the cell: along directions
# where the tessellation is periodic but the cell is not, its translates clipped to the cell
function _cell_pieces(t::Tessellation{D}, P) where {D}
    cut = [d for d in 1:D if t.periodic[d] && !t.cell.periodic[d]]
    isempty(cut) && return [P]
    L = t.cell.size
    eps = 1e-11 * maximum(L)
    out = Any[]
    for s in Iterators.product(ntuple(d -> d in cut ? (-1:1) : (0:0), D)...)
        Q = D == 3 ? _Poly3([[p .+ s .* L for p in f] for f in P.faces], copy(P.labels)) :
                     _Poly2([p .+ s .* L for p in P.pts], copy(P.labels))
        for d in cut, (n, c) in ((ntuple(k -> k == d ? 1.0 : 0.0, D), L[d]), (ntuple(k -> k == d ? -1.0 : 0.0, D), 0.0))
            _clip!(Q, n, c, (0, ntuple(_ -> 0, D)), eps)
        end
        (D == 3 ? length(Q.faces) >= 4 : length(Q.pts) >= 3) || continue
        _poly_volume(Q) > 1e-12 * prod(L) && push!(out, Q)
    end
    return out
end

"""
    pores(t; wall) -> Vector{Polyhedron} (3D) or Vector{Polygon} (2D)

The grains shrunk by `wall / 2` on every face: the pores of a closed-cell foam whose walls
(the rest of the cell) have thickness `wall`. Grains that vanish are skipped with a warning.
Use `RVE(cell, Phase(:pores, pores(t; wall); void = true); matrix = :walls)` for empty pores,
or without `void` for a two-phase material. [`wall_thickness`](@ref) gives the wall for a
relative density. To refine the pore surfaces, make the walls an explicit phase:
`[Phase(:walls, Box(cell.origin, cell.size); refine = w / 2), Phase(:pores, pores(t; wall = w); void = true)]`
with `matrix = nothing`. At free faces of the cell the wall is `wall / 2` thick.
"""
function pores(t::Tessellation{D}; wall) where {D}
    w = Float64(wall)
    w >= 0 || throw(ArgumentError("pores: wall must be >= 0"))
    out = D == 3 ? Polyhedron[] : Polygon[]
    lost = 0
    V = prod(t.cell.size)
    for i in eachindex(t.loops)
        P = _shrunk(t, i, w / 2)
        if P === nothing || _poly_volume(P) <= 1e-9 * V
            lost += 1                   # vanished (or too small to be a valid body)
            continue
        end
        for Q in _cell_pieces(t, P)
            try
                if D == 3
                    v, f = _indexed(Q, 1e-9 * maximum(t.cell.size))
                    push!(out, Polyhedron([t.cell.origin .+ p for p in v], f))
                else
                    push!(out, Polygon([t.cell.origin .+ p for p in Q.pts]))
                end
            catch e
                e isa ArgumentError || rethrow()
                lost += 1               # a sliver of a nearly vanished pore
            end
        end
    end
    lost > 0 && @warn "pores: $lost grains vanish with wall = $w"
    return out
end

"""
    wall_thickness(t, density) -> wall

Wall thickness of the closed-cell foam [`pores`](@ref)`(t; wall)` with the given relative
density (volume fraction of the walls), exact for the tessellation (bisection on the pore
volumes).
"""
function wall_thickness(t::Tessellation{D}, density::Real) where {D}
    0 < density < 1 || throw(ArgumentError("wall_thickness: density must be in (0, 1)"))
    V = prod(t.cell.size)
    pv(P) = P === nothing ? 0.0 : sum(_poly_volume, _cell_pieces(t, P); init=0.0)
    ρ(w) = 1 - sum(pv(_shrunk(t, i, w / 2)) for i in eachindex(t.loops)) / V
    lo, hi = 0.0, minimum(t.cell.size)
    for _ in 1:200
        m = (lo + hi) / 2
        ρ(m) < density ? (lo = m) : (hi = m)
        hi - lo < 1e-13 * maximum(t.cell.size) && break
    end
    return (lo + hi) / 2
end
