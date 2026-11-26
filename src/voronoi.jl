# Periodic Voronoi tessellations (pure Julia): every grain is the box seed ± L/2 clipped by the
# bisectors of the other seeds and their periodic images; the grains are then welded into one
# periodic cell complex so that neighbouring grains share bit-identical vertices.

const _VQ = 1e-8                     # weld tolerance (fractional); merged vertices are < 1e-7 apart, as Polyhedron requires

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

function _voronoi_cells(X::Vector{NTuple{D,Float64}}, L::NTuple{D,Float64}, per::NTuple{D,Bool}) where {D}
    N = length(X)
    eps = 1e-11 * maximum(L)
    cells = Any[]
    for i in 1:N
        P = D == 3 ? _box3(X[i], L, per, i) : _box2(X[i], L, per, i)
        R0 = maximum(p -> _norm(p .- X[i]), _vertices(P))
        K = ntuple(d -> per[d] ? ceil(Int, 2R0 / L[d]) : 0, D)
        cand = Tuple{Float64,Int,NTuple{D,Int}}[]
        for j in 1:N, k in Iterators.product(ntuple(d -> -K[d]:K[d], D)...)
            j == i && all(==(0), k) && continue
            y = X[j] .+ k .* L
            r = _norm(y .- X[i])
            r < 2R0 + eps && push!(cand, (r, j, k))
        end
        sort!(cand)
        R = R0
        for (r, j, k) in cand
            r > 2R + eps && break
            y = X[j] .+ k .* L
            n = (y .- X[i]) ./ r
            if _clip!(P, n, _dot(n, (X[i] .+ y) ./ 2), (j, k), eps)
                R = maximum(p -> _norm(p .- X[i]), _vertices(P))
            end
        end
        push!(cells, P)
    end
    return cells
end

# ---------------------------------------------------------------- welded complex

"""
    Tessellation{D}

Periodic Voronoi tessellation of a cell (see [`voronoi`](@ref)): seeds, welded vertices
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
end

Base.show(io::IO, t::Tessellation{D}) where {D} =
    print(io, "Tessellation{$D}($(length(t.seeds)) grains, $(length(t.vertices)) vertices)")

function _weld(cells, L::NTuple{D,Float64}, per) where {D}
    V = Vec{D}[]
    bins = Dict{NTuple{D,Int},Vector{Int}}()
    nb = round(Int, 1 / _VQ)
    function vid(p)
        f = p ./ L
        w = ntuple(d -> per[d] ? f[d] - floor(f[d]) : f[d], D)
        key = ntuple(d -> round(Int, w[d] / _VQ), D)
        for s in Iterators.product(ntuple(_ -> -1:1, D)...)
            kk = ntuple(d -> per[d] ? mod(key[d] + s[d], nb) : key[d] + s[d], D)
            for id in get(bins, kk, Int[])
                δ = ntuple(D) do d
                    x = f[d] - V[id][d]
                    per[d] ? x - round(x) : x
                end
                if maximum(abs, δ) <= 10_VQ
                    return id, ntuple(d -> per[d] ? round(Int, f[d] - V[id][d]) : 0, D)
                end
            end
        end
        push!(V, w)
        push!(get!(bins, ntuple(d -> per[d] ? mod(key[d], nb) : key[d], D), Int[]), length(V))
        return length(V), ntuple(d -> per[d] ? round(Int, f[d] - w[d]) : 0, D)
    end
    loops = Vector{Vector{Vector{Tuple{Int,NTuple{D,Int}}}}}()
    nbrs = Vector{Vector{Tuple{Int,NTuple{D,Int}}}}()
    # loops of welded vertices; consecutive repeats (edges below the weld tolerance) removed
    function wloop(pts)
        l = [vid(p) for p in pts]
        keep = [l[k] != l[mod1(k + 1, length(l))] for k in eachindex(l)]
        return l[keep], keep
    end
    for P in cells
        if D == 3
            F = Vector{Tuple{Int,NTuple{D,Int}}}[]; lab = Tuple{Int,NTuple{D,Int}}[]
            for (f, l) in zip(P.faces, P.labels)
                w, _ = wloop(f)
                length(w) >= 3 && (push!(F, w); push!(lab, l))
            end
            push!(loops, F); push!(nbrs, lab)
        else
            w, keep = wloop(P.pts)
            push!(loops, [w]); push!(nbrs, P.labels[keep])
        end
    end
    return V, loops, nbrs
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
block cut from a periodic tessellation. `rng`: an `AbstractRNG` or an integer seed.

Then `RVE(cell, grain_phases(t); matrix = nothing)` gives a polycrystal with one phase per grain
and the grain boundaries as interface groups.
"""
function voronoi(cell::Cell{D}; n=nothing, seeds=nothing, min_distance=:auto, lloyd=0,
                 periodic=cell.periodic, rng=Random.default_rng(), max_tries=50_000) where {D}
    (n === nothing) != (seeds === nothing) || throw(ArgumentError("voronoi: give exactly one of n or seeds"))
    per = periodic isa Bool ? ntuple(_ -> periodic, D) : NTuple{D,Bool}(periodic)
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
        n isa Integer && n >= 1 || throw(ArgumentError("voronoi: n must be a positive integer"))
        V = prod(L)
        dmin = min_distance === :auto ? 0.5 * (V / n)^(1 / D) :
               min_distance isa Real && min_distance >= 0 ? Float64(min_distance) :
               throw(ArgumentError("voronoi: min_distance must be :auto or a number >= 0"))
        r = _rng(rng)
        X = if dmin == 0
            [ntuple(d -> rand(r) * L[d], D) for _ in 1:n]
        else
            2 * (dmin / 2) < minimum(L) || throw(ArgumentError("voronoi: min_distance must be smaller than the cell"))
            Y = _rsa(r, L, per, fill(dmin / 2, n), 0.0, 0.0, false, max_tries)
            length(Y) == n || throw(ArgumentError(@sprintf("voronoi: could not place %d seeds %.3g apart (placed %d); reduce min_distance",
                                                           n, dmin, length(Y))))
            Y
        end
    end
    for i in eachindex(X), j in i+1:length(X)
        δ = ntuple(d -> (x = X[i][d] - X[j][d]; per[d] ? x - L[d] * round(x / L[d]) : x), D)
        _norm(δ) > _VQ * maximum(L) || throw(ArgumentError("voronoi: seeds $i and $j coincide"))
    end
    t = _tessellate(cell, X, L, per)
    for _ in 1:lloyd
        c = [_grain_centroid(t, i) for i in eachindex(X)]                   # cell-local
        X = [ntuple(d -> per[d] ? c[i][d] - L[d] * floor(c[i][d] / L[d]) : c[i][d], D) for i in eachindex(c)]
        t = _tessellate(cell, X, L, per)
    end
    return t
end

function _tessellate(cell::Cell{D}, X, L, per) where {D}
    cells = _voronoi_cells(X, L, per)
    V, loops, nbrs = _weld(cells, L, per)
    t = Tessellation{D}(cell, per, [x ./ L for x in X], V, loops, nbrs)
    _check_faces(t)
    return t
end

# cell-local coordinates of a welded vertex with an offset
_vpos(t::Tessellation{D}, (v, o)) where {D} = ntuple(d -> (t.vertices[v][d] + o[d]) * t.cell.size[d], D)

function _grain_geometry(t::Tessellation{3}, i)
    ids = Dict{Tuple{Int,NTuple{3,Int}},Int}()
    verts = Vec{3}[]
    faces = Vector{Int}[]
    for f in t.loops[i]
        push!(faces, [get!(ids, vo) do
                          push!(verts, _vpos(t, vo)); length(verts)
                      end for vo in f])
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
grains(t::Tessellation{3}) = [let (v, f) = _grain_geometry(t, i)
                                  Polyhedron([t.cell.origin .+ p for p in v], f)
                              end for i in eachindex(t.loops)]
grains(t::Tessellation{2}) = [Polygon([t.cell.origin .+ _vpos(t, vo) for vo in only(t.loops[i])]) for i in eachindex(t.loops)]

"""
    grain_volumes(t) -> Vector{Float64}

Volume (area in 2D) of each grain; they add up to the volume of the cell.
"""
grain_volumes(t::Tessellation{3}) = [abs(_polyhedron_volume(_grain_geometry(t, i)...)) for i in eachindex(t.loops)]
grain_volumes(t::Tessellation{2}) = [abs(_polygon_area([_vpos(t, vo) for vo in only(t.loops[i])])) for i in eachindex(t.loops)]

"""
    grain_phases(t; prefix = "grain", groups = nothing) -> Vector{Phase}

One phase per grain, named `grain_1`, `grain_2`, …; or, with `groups` (a vector of vectors of
grain indices, or a vector with the group number of each grain), one fused phase per group,
named `grain_1`, … in the order of the groups. Use with `matrix = nothing`:
`RVE(cell, grain_phases(t); matrix = nothing)`; the grain boundaries are the interface groups
`interface_grain_i_grain_j`.
"""
function grain_phases(t::Tessellation{D}; prefix="grain", groups=nothing) where {D}
    g = grains(t)
    groups === nothing && return [Phase(Symbol(prefix, "_", i), g[i]) for i in eachindex(g)]
    sets = if eltype(groups) <: Integer
        length(groups) == length(g) || throw(ArgumentError("grain_phases: one group number per grain"))
        [findall(==(k), groups) for k in sort(unique(groups))]
    else
        [collect(Int, s) for s in groups]
    end
    sort(reduce(vcat, sets)) == collect(eachindex(g)) ||
        throw(ArgumentError("grain_phases: every grain must be in exactly one group"))
    return [Phase(Symbol(prefix, "_", k), g[s]; fuse=length(s) > 1) for (k, s) in enumerate(sets)]
end
