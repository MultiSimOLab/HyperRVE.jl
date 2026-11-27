# Periodic Voronoi tessellations (pure Julia): every grain is the box seed ± L/2 clipped by the
# bisectors of the other seeds and their periodic images; the grains are then welded into one
# periodic cell complex so that neighbouring grains share bit-identical vertices.

const _VQ = 2e-7                     # weld tolerance (fractional): vertices left distinct are farther apart than Polyhedron's 1e-7

# ---------------------------------------------------------------- 3D clipping

# convex polyhedron: faces as point loops, each with a label (neighbour j, lattice offset k);
# j = 0 is a wall of a non-periodic direction
mutable struct _Poly3
    faces::Vector{Vector{NTuple{3,Float64}}}
    labels::Vector{Tuple{Int,NTuple{3,Int}}}
end

function _box3(x, L, per, i)
    lo = ntuple(d -> per[d] ? x[d] - L[d] / 2 : 0.0, 3)
    hi = ntuple(d -> per[d] ? x[d] + L[d] / 2 : L[d], 3)
    c(a, b, cc) = (a == 0 ? lo[1] : hi[1], b == 0 ? lo[2] : hi[2], cc == 0 ? lo[3] : hi[3])
    F = [[c(0, 0, 0), c(0, 0, 1), c(0, 1, 1), c(0, 1, 0)], [c(1, 0, 0), c(1, 1, 0), c(1, 1, 1), c(1, 0, 1)],
         [c(0, 0, 0), c(1, 0, 0), c(1, 0, 1), c(0, 0, 1)], [c(0, 1, 0), c(0, 1, 1), c(1, 1, 1), c(1, 1, 0)],
         [c(0, 0, 0), c(0, 1, 0), c(1, 1, 0), c(1, 0, 0)], [c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)]]
    lab = Tuple{Int,NTuple{3,Int}}[]
    for d in 1:3, s in (-1, 1)
        e = ntuple(k -> k == d ? s : 0, 3)
        push!(lab, (per[d] ? i : 0, e))
    end
    return _Poly3(F, lab)
end

_dot(a, b) = sum(a .* b)
_cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

function _push_unique!(v, p, tol)
    any(q -> _norm(q .- p) <= tol, v) || push!(v, p)
    return v
end

"Clip by the half-space n·x ≤ c; returns whether anything was removed."
function _clip!(P::_Poly3, n, c, label, eps)
    any(f -> any(p -> _dot(n, p) - c > eps, f), P.faces) || return false
    faces = Vector{NTuple{3,Float64}}[]
    labels = Tuple{Int,NTuple{3,Int}}[]
    cap = NTuple{3,Float64}[]
    for (f, l) in zip(P.faces, P.labels)
        out = NTuple{3,Float64}[]
        m = length(f)
        for k in 1:m
            p, q = f[k], f[mod1(k + 1, m)]
            dp, dq = _dot(n, p) - c, _dot(n, q) - c
            if dp <= eps
                push!(out, p)
                abs(dp) <= eps && _push_unique!(cap, p, 10eps)
            end
            if (dp < -eps && dq > eps) || (dp > eps && dq < -eps)
                t = dp / (dp - dq)
                r = p .+ t .* (q .- p)
                push!(out, r)
                _push_unique!(cap, r, 10eps)
            end
        end
        g = NTuple{3,Float64}[]
        for p in out
            _push_unique!(g, p, 10eps)
        end
        if length(g) >= 3
            push!(faces, g); push!(labels, l)
        end
    end
    if length(cap) >= 3
        ctr = reduce(.+, cap) ./ length(cap)
        u = _cross(n, abs(n[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0))
        u = u ./ _norm(u)
        w = _cross(n, u)
        sort!(cap; by=p -> atan(_dot(p .- ctr, w), _dot(p .- ctr, u)))
        push!(faces, cap); push!(labels, label)
    end
    P.faces, P.labels = faces, labels
    return true
end

# ---------------------------------------------------------------- 2D clipping

# convex polygon: loop of points, label of the edge from point k to point k + 1
mutable struct _Poly2
    pts::Vector{NTuple{2,Float64}}
    labels::Vector{Tuple{Int,NTuple{2,Int}}}
end

function _box2(x, L, per, i)
    lo = ntuple(d -> per[d] ? x[d] - L[d] / 2 : 0.0, 2)
    hi = ntuple(d -> per[d] ? x[d] + L[d] / 2 : L[d], 2)
    pts = [(lo[1], lo[2]), (hi[1], lo[2]), (hi[1], hi[2]), (lo[1], hi[2])]
    lab = [(per[2] ? i : 0, (0, -1)), (per[1] ? i : 0, (1, 0)), (per[2] ? i : 0, (0, 1)), (per[1] ? i : 0, (-1, 0))]
    return _Poly2(pts, lab)
end

function _clip!(P::_Poly2, n, c, label, eps)
    any(p -> _dot(n, p) - c > eps, P.pts) || return false
    # Sutherland–Hodgman; the label of a point is the label of the edge that starts at it
    pts = NTuple{2,Float64}[]; labs = Tuple{Int,NTuple{2,Int}}[]
    m = length(P.pts)
    for k in 1:m
        p, q = P.pts[k], P.pts[mod1(k + 1, m)]
        dp, dq = _dot(n, p) - c, _dot(n, q) - c
        if dp <= eps
            push!(pts, p); push!(labs, P.labels[k])
            if dq > eps                                   # leaves: the cut edge starts at the exit point
                t = dp / (dp - dq)
                push!(pts, p .+ t .* (q .- p)); push!(labs, label)
            end
        elseif dq <= eps                                  # enters: the rest of edge k starts at the entry point
            t = dp / (dp - dq)
            push!(pts, p .+ t .* (q .- p)); push!(labs, P.labels[k])
        end
    end
    keep = [_norm(pts[k] .- pts[mod1(k + 1, length(pts))]) > 10eps for k in eachindex(pts)]
    P.pts, P.labels = pts[keep], labs[keep]
    return true
end

_vertices(P::_Poly3) = unique(reduce(vcat, P.faces))
_vertices(P::_Poly2) = P.pts

# ---------------------------------------------------------------- cells of all seeds

# Power (Laguerre) cells: the bisector of seeds i and j (weights w) is the plane at distance
# (r² + w_i − w_j) / (2r) from x_i along the direction to x_j; w = 0 gives Voronoi cells. The
# candidates are sorted by that distance and the loop stops when the plane lies beyond the
# farthest vertex of the cell. A cell can be empty (a seed dominated by its neighbours).
function _voronoi_cells(X::Vector{NTuple{D,Float64}}, L::NTuple{D,Float64}, per::NTuple{D,Bool}, epsrel=1e-11,
                        w::AbstractVector{Float64}=zeros(length(X))) where {D}
    N = length(X)
    eps = epsrel * maximum(L)
    wmax = maximum(w)
    cells = Any[]
    for i in 1:N
        P = D == 3 ? _box3(X[i], L, per, i) : _box2(X[i], L, per, i)
        R0 = maximum(p -> _norm(p .- X[i]), _vertices(P))
        # a plane can cut the box only if (r² + w_i − w_j)/(2r) <= R0
        rmax = R0 + sqrt(R0^2 + max(wmax - w[i], 0.0))
        cand = Tuple{Float64,Float64,Int,NTuple{D,Int}}[]
        for j in 1:N
            kr = ntuple(d -> per[d] ? (ceil(Int, (X[i][d] - rmax - X[j][d]) / L[d]):floor(Int, (X[i][d] + rmax - X[j][d]) / L[d])) : (0:0), D)
            for k in Iterators.product(kr...)
                j == i && all(==(0), k) && continue
                y = X[j] .+ k .* L
                r = _norm(y .- X[i])
                r < rmax + eps && push!(cand, ((r^2 + w[i] - w[j]) / (2r), r, j, k))
            end
        end
        sort!(cand)
        R = R0
        for (dpl, r, j, k) in cand
            dpl > R + eps && break
            y = X[j] .+ k .* L
            n = (y .- X[i]) ./ r
            if _clip!(P, n, _dot(n, X[i]) + dpl, (j, k), eps)
                (D == 3 ? isempty(P.faces) : length(P.pts) < 3) && throw(_EmptyCell(i))
                R = maximum(p -> _norm(p .- X[i]), _vertices(P))
            end
        end
        push!(cells, P)
    end
    return cells
end

struct _EmptyCell <: Exception
    i::Int
end

# ---------------------------------------------------------------- welded complex

"""
    Tessellation{D}

Periodic Voronoi or Laguerre tessellation of a cell (see [`voronoi`](@ref) and
[`laguerre`](@ref)): seeds, power weights (zero for Voronoi), welded vertices
(fractional coordinates in `[0, 1)` along the wrapped directions) and, for each grain, its
faces (3D: loops of `(vertex, offset)`; 2D: one loop) with the neighbour across each face.
Use [`grains`](@ref), [`grain_phases`](@ref) and [`grain_volumes`](@ref).
"""
struct Tessellation{D}
    cell::Cell{D}
    periodic::NTuple{D,Bool}
    seeds::Vector{Vec{D}}
    vertices::Vector{Vec{D}}
    loops::Vector{Vector{Vector{Tuple{Int,NTuple{D,Int}}}}}     # grain → faces (2D: one loop) → (vertex, offset)
    neighbours::Vector{Vector{Tuple{Int,NTuple{D,Int}}}}        # grain → face → (neighbour, offset); 0 = wall
    moved::Vector{Bool}                                         # vertex displaced by the weld (its faces are triangulated)
    weights::Vector{Float64}                                    # power weights (zeros: Voronoi)
end

Base.show(io::IO, t::Tessellation{D}) where {D} =
    print(io, "Tessellation{$D}($(length(t.seeds)) grains, $(length(t.vertices)) vertices)")

function _weld(cells, L::NTuple{D,Float64}, per, VQ=_VQ) where {D}
    # all raw points of all grains, fractional
    raw = NTuple{D,Float64}[]
    for P in cells
        for p in (D == 3 ? Iterators.flatten(P.faces) : P.pts)
            push!(raw, p ./ L)
        end
    end
    wrapf(f) = ntuple(d -> per[d] ? f[d] - floor(f[d]) : f[d], D)
    pdist(a, b) = maximum(ntuple(D) do d
        x = a[d] - b[d]
        abs(per[d] ? x - round(x) : x)
    end)
    # union–find of the points closer than VQ (transitive, the same for every grain)
    parent = collect(eachindex(raw))
    findr(i) = (while parent[i] != i; parent[i] = parent[parent[i]]; i = parent[i]; end; i)
    nb = max(1, floor(Int, 1 / VQ))
    key(w) = ntuple(d -> floor(Int, w[d] / VQ), D)
    bins = Dict{NTuple{D,Int},Vector{Int}}()
    W = [wrapf(f) for f in raw]
    for (i, w) in enumerate(W)
        k = key(w)
        for s in Iterators.product(ntuple(_ -> -1:1, D)...)
            kk = ntuple(d -> per[d] ? mod(k[d] + s[d], nb) : k[d] + s[d], D)
            for j in get(bins, kk, Int[])
                if pdist(w, W[j]) <= VQ
                    a, b = findr(i), findr(j)
                    a != b && (parent[max(a, b)] = min(a, b))
                end
            end
        end
        push!(get!(bins, ntuple(d -> per[d] ? mod(k[d], nb) : k[d], D), Int[]), i)
    end
    # one vertex per class at the position of its root; offsets of every raw point
    vid = Dict{Int,Int}()
    V = Vec{D}[]
    spread = Float64[]
    for i in eachindex(raw)
        r = findr(i)
        id = get!(vid, r) do
            push!(V, W[r]); push!(spread, 0.0); length(V)
        end
        spread[id] = max(spread[id], pdist(W[i], W[r]))
    end
    moved = spread .> 1e-12                     # vertices displaced by the weld
    off(i) = ntuple(d -> per[d] ? round(Int, raw[i][d] - V[vid[findr(i)]][d]) : 0, D)
    loops = Vector{Vector{Vector{Tuple{Int,NTuple{D,Int}}}}}()
    nbrs = Vector{Vector{Tuple{Int,NTuple{D,Int}}}}()
    n = 0
    function wloop(m)
        l = [(vid[findr(n + k)], off(n + k)) for k in 1:m]
        n += m
        keep = [l[k] != l[mod1(k + 1, length(l))] for k in eachindex(l)]
        return l[keep], keep
    end
    for P in cells
        if D == 3
            F = Vector{Tuple{Int,NTuple{D,Int}}}[]; lab = Tuple{Int,NTuple{D,Int}}[]
            for (f, l) in zip(P.faces, P.labels)
                w, _ = wloop(length(f))
                length(unique(w)) >= 3 && (push!(F, w); push!(lab, l))
            end
            push!(loops, F); push!(nbrs, lab)
        else
            w, keep = wloop(length(P.pts))
            push!(loops, [w]); push!(nbrs, P.labels[keep])
        end
    end
    return V, loops, nbrs, moved
end

# every face (a → (b, k)) must appear as (b → (a, −k)) with the same vertices shifted by −k
function _check_faces(t::Tessellation{D}) where {D}
    D == 3 || return nothing
    faces = Dict{Tuple{Int,Int,NTuple{3,Int}},Set{Tuple{Int,NTuple{3,Int}}}}()
    for a in eachindex(t.loops), (f, (b, k)) in zip(t.loops[a], t.neighbours[a])
        b == 0 && continue
        faces[(a, b, k)] = Set(f)
    end
    for ((a, b, k), s) in faces
        o = get(faces, (b, a, .-k), nothing)
        o === nothing && error("voronoi: internal error, the face between grains $a and $b has one side only; " *
                               "try other seeds (degenerate configuration)")
        Set((v, off .- k) for (v, off) in s) == o ||
            error("voronoi: internal error, grains $a and $b do not share the same face vertices")
    end
    return nothing
end

"""
    voronoi(cell; n = nothing, seeds = nothing, min_distance = :auto, lloyd = 0,
            periodic = cell.periodic, rng = Random.default_rng(), max_tries = 50_000) -> Tessellation

Periodic Voronoi tessellation of the cell with `n` random seeds or the given `seeds`
(physical points). Random seeds are at least `min_distance` apart (`:auto` = half the mean
grain size, `0` = Poisson–Voronoi). `lloyd = k` moves the seeds `k` times to the centroids of
their grains (more uniform grain sizes). `periodic = true` on a cell with free faces gives a
block cut from a periodic tessellation (its grains are clipped at the free faces).
`rng`: an `AbstractRNG` or an integer seed; `max_tries`: attempts of the seed placement.

Vertices closer than 2e-7 of the cell are merged (the faces they touch are then split into
triangles, identically in both grains). Seeds on a slightly perturbed regular lattice give
grain faces nearly coplanar with the cell faces, which Gmsh may fail to make periodic.

Then `RVE(cell, grain_phases(t); matrix = nothing)` gives a polycrystal with one phase per grain
and the grain boundaries as interface groups.
"""
function voronoi(cell::Cell{D}; n=nothing, seeds=nothing, min_distance=:auto, lloyd=0,
                 periodic=cell.periodic, rng=Random.default_rng(), max_tries=50_000) where {D}
    (n === nothing) != (seeds === nothing) || throw(ArgumentError("voronoi: give exactly one of n or seeds"))
    per = periodic isa Bool ? ntuple(_ -> periodic, D) :
          (length(periodic) == D && all(x -> x isa Bool, periodic)) ? NTuple{D,Bool}(Tuple(periodic)) :
          throw(ArgumentError("voronoi: periodic must be a Bool or $D Bools"))
    L = ntuple(d -> Float64(cell.size[d]), D)
    lloyd isa Integer && lloyd >= 0 || throw(ArgumentError("voronoi: lloyd must be an integer >= 0"))
    if seeds !== nothing
        X = NTuple{D,Float64}[]
        for s in seeds
            p = _vec(Val(D), s, "voronoi seed") .- cell.origin
            p = ntuple(D) do d
                per[d] ? p[d] - L[d] * floor(p[d] / L[d]) :
                (0 <= p[d] <= L[d] ? p[d] : throw(ArgumentError("voronoi: seed $s is outside the cell along the non-periodic direction $d")))
            end
            push!(X, p)
        end
        isempty(X) && throw(ArgumentError("voronoi: no seeds"))
    else
        n isa Integer && !(n isa Bool) && n >= 1 || throw(ArgumentError("voronoi: n must be a positive integer"))
        V = prod(L)
        # :auto: half the mean grain size, but never more than fits in the cell
        dmin = min_distance === :auto ? (n == 1 ? 0.0 : min(0.5 * (V / n)^(1 / D), 0.45 * minimum(L))) :
               min_distance isa Real && min_distance >= 0 ? Float64(min_distance) :
               throw(ArgumentError("voronoi: min_distance must be :auto or a number >= 0"))
        r = _rng(rng)
        X = if dmin == 0
            [ntuple(d -> rand(r) * L[d], D) for _ in 1:n]
        else
            dmin < minimum(L) || throw(ArgumentError("voronoi: min_distance must be smaller than the cell"))
            Y = _rsa(r, L, per, fill(dmin / 2, n), 0.0, 0.0, false, max_tries)
            length(Y) == n || throw(ArgumentError(@sprintf("voronoi: could not place %d seeds %.3g apart (placed %d); reduce min_distance (about %.3g is feasible)",
                                                           n, dmin, length(Y), (_SAT[D] * V / (n * _ballvol(Val(D), 1.0)))^(1 / D) * 2 * 0.9)))
            Y
        end
    end
    for i in eachindex(X), j in i+1:length(X)
        δ = ntuple(d -> (x = X[i][d] - X[j][d]; per[d] ? x - L[d] * round(x / L[d]) : x), D)
        _norm(δ) > 1e-6 * maximum(L) || throw(ArgumentError("voronoi: seeds $i and $j coincide (closer than 1e-6 of the cell)"))
    end
    t = _tessellate(cell, X, L, per)
    for _ in 1:lloyd
        c = [_grain_centroid(t, i) for i in eachindex(X)]                   # cell-local
        X = [ntuple(d -> per[d] ? c[i][d] - L[d] * floor(c[i][d] / L[d]) : c[i][d], D) for i in eachindex(c)]
        t = _tessellate(cell, X, L, per)
    end
    return t
end

# Nearly degenerate seeds (a slightly perturbed regular lattice) make the in/on/out
# classification of nearly coincident planes differ between neighbouring grains; a larger
# classification tolerance treats them consistently as degenerate.
function _tessellate(cell::Cell{D}, X, L, per, w=zeros(length(X)); what="voronoi") where {D}
    msg = ""
    for (epsrel, vq) in ((1e-11, _VQ), (1e-9, _VQ), (1e-7, 1e-6), (1e-5, 1e-4))
        cells = try
            _voronoi_cells(X, L, per, epsrel, w)
        catch e
            e isa _EmptyCell || rethrow()
            throw(ArgumentError("$what: the cell of seed $(e.i) is empty (its weight is too small for its neighbours)"))
        end
        V, loops, nbrs, moved = _weld(cells, L, per, vq)
        t = Tessellation{D}(cell, per, [x ./ L for x in X], V, loops, nbrs, moved, collect(Float64, w))
        try
            _check_faces(t)
            grains(t)                                   # every grain must be a valid Polyhedron / Polygon
            return t
        catch e
            e isa Union{ErrorException,ArgumentError} || rethrow()
            msg = sprint(showerror, e)
        end
    end
    throw(ArgumentError("$what: the seeds are too close to a degenerate configuration (e.g. a regular lattice " *
                        "perturbed by much less than its spacing); use the exact lattice or perturb it by more than " *
                        "1e-3 of the spacing ($msg)"))
end

# cell-local coordinates of a welded vertex with an offset
_vpos(t::Tessellation{D}, (v, o)) where {D} = ntuple(d -> (t.vertices[v][d] + o[d]) * t.cell.size[d], D)

# a face with a vertex moved by the weld is no longer planar: it is split into a fan of
# triangles from its smallest (vertex, offset), which is the same vertex seen from both grains
# (their offsets differ by one lattice shift), so the two sides get the same triangles
function _face_pieces(t::Tessellation, f)
    any(vo -> t.moved[vo[1]], f) && length(f) > 3 || return [f]
    a = argmin(f)
    g = circshift(f, 1 - a)
    return [[g[1], g[k], g[k+1]] for k in 2:length(g)-1]
end

function _grain_geometry(t::Tessellation{3}, i)
    ids = Dict{Tuple{Int,NTuple{3,Int}},Int}()
    verts = Vec{3}[]
    faces = Vector{Int}[]
    for f in t.loops[i], piece in _face_pieces(t, f)
        push!(faces, [get!(ids, vo) do
                          push!(verts, _vpos(t, vo)); length(verts)
                      end for vo in piece])
    end
    return verts, faces
end

function _grain_centroid(t::Tessellation{3}, i)
    v, f = _grain_geometry(t, i)
    return _polyhedron_centroid(v, f)
end
_grain_centroid(t::Tessellation{2}, i) = _polygon_centroid([_vpos(t, vo) for vo in only(t.loops[i])])

"""
    grains(t::Tessellation) -> Vector{Polyhedron} (3D) or Vector{Polygon} (2D)

The whole grains in physical coordinates (a grain may stick out of the cell; its periodic
images complete the cell when it is meshed).
"""
grains(t::Tessellation) = [grains(t, i) for i in eachindex(t.loops)]
grains(t::Tessellation{3}, i::Integer) = let (v, f) = _grain_geometry(t, i)
    Polyhedron([t.cell.origin .+ p for p in v], f)
end
grains(t::Tessellation{2}, i::Integer) = Polygon([t.cell.origin .+ _vpos(t, vo) for vo in only(t.loops[i])])

"""
    grain_volumes(t) -> Vector{Float64}

Volume (area in 2D) of each grain; they add up to the volume of the cell.
"""
grain_volumes(t::Tessellation{3}) = [abs(_polyhedron_volume(_grain_geometry(t, i)...)) for i in eachindex(t.loops)]
grain_volumes(t::Tessellation{2}) = [abs(_polygon_area([_vpos(t, vo) for vo in only(t.loops[i])])) for i in eachindex(t.loops)]

# polyhedron / polygon of a clipping body (cell-local) in physical coordinates
function _indexed(P::_Poly3, tol)
    verts = NTuple{3,Float64}[]; faces = Vector{Int}[]
    for f in P.faces
        push!(faces, [begin
                          k = findfirst(q -> _norm(q .- p) <= tol, verts)
                          k === nothing ? (push!(verts, p); length(verts)) : k
                      end for p in f])
    end
    return verts, faces
end

# Items of grain i for a phase: the whole grain, or, along directions where the tessellation
# is periodic but the cell is not (a block cut from a periodic tessellation), its translates
# clipped to the cell
function _grain_items(t::Tessellation{D}, i) where {D}
    cut = [d for d in 1:D if t.periodic[d] && !t.cell.periodic[d]]
    isempty(cut) && return Geometry{D}[grains(t, i)]
    L = t.cell.size
    eps = 1e-11 * maximum(L)
    out = Geometry{D}[]
    for s in Iterators.product(ntuple(d -> d in cut ? (-1:1) : (0:0), D)...)
        P = if D == 3
            v, f = _grain_geometry(t, i)
            _Poly3([[v[k] .+ s .* L for k in face] for face in f], [(0, (0, 0, 0)) for _ in f])
        else
            _Poly2([_vpos(t, vo) .+ s .* L for vo in only(t.loops[i])], [(0, (0, 0)) for _ in only(t.loops[i])])
        end
        for d in cut, (n, c) in ((ntuple(k -> k == d ? 1.0 : 0.0, D), L[d]), (ntuple(k -> k == d ? -1.0 : 0.0, D), 0.0))
            _clip!(P, n, c, (0, ntuple(_ -> 0, D)), eps)
        end
        if D == 3
            length(P.faces) >= 4 || continue
            v, f = _indexed(P, 1e-9 * maximum(L))
            abs(_polyhedron_volume(v, f)) > 1e-12 * prod(L) || continue
            push!(out, Polyhedron([t.cell.origin .+ p for p in v], f))
        else
            length(P.pts) >= 3 && abs(_polygon_area(P.pts)) > 1e-12 * prod(L) || continue
            push!(out, Polygon([t.cell.origin .+ p for p in P.pts]))
        end
    end
    return out
end

"""
    grain_phases(t; prefix = "grain", groups = nothing, fuse = groups !== nothing, h = nothing, refine = nothing)

One phase per grain, named `grain_1`, `grain_2`, …; or, with `groups` (one label per grain),
one phase per distinct label: `<prefix>_<label>` for integer labels, the label itself for
symbols and strings. With `fuse = true` (the default with groups) the grains of a phase are
united, so no surfaces remain between them. `h` and `refine` are passed to every phase (e.g.
`refine = 0.01` refines along the grain boundaries). Use without matrix:
`RVE(cell, grain_phases(t); matrix = nothing)`; the grain boundaries are the interface groups
`interface_grain_i_grain_j`. A tessellation made with `periodic = true` on a cell with free
faces is cut at those faces (its grains are clipped there).
"""
function grain_phases(t::Tessellation{D}; prefix="grain", groups=nothing, fuse=groups !== nothing,
                      h=nothing, refine=nothing) where {D}
    N = length(t.loops)
    items(i) = _grain_items(t, i)
    groups === nothing && return [Phase(Symbol(prefix, "_", i), items(i); h, refine, fuse=false) for i in 1:N]
    (groups isa AbstractVector && length(groups) == N) ||
        throw(ArgumentError("grain_phases: groups must give one label per grain ($N)"))
    all(x -> (x isa Integer && !(x isa Bool)) || x isa Symbol || x isa AbstractString, groups) ||
        throw(ArgumentError("grain_phases: group labels must be integers, symbols or strings"))
    labels = unique(groups)
    all(x -> x isa Integer, labels) && sort!(labels)
    name(l) = l isa Integer ? Symbol(prefix, "_", l) : Symbol(l)
    return [Phase(name(l), reduce(vcat, [items(i) for i in 1:N if groups[i] == l]); h, refine, fuse) for l in labels]
end

