# Rounded pores of closed-cell foams: each pore is its shrunk grain intersected with a ball
# centred at the grain centroid, of radius (equivalent radius of the shrunk grain) / rounding.
# Volumes of a convex body ∩ ball are integrated along the rays from the centre.

# quadrature directions: spherical Fibonacci points (3D), equally spaced angles (2D)
const _RAYS3 = let M = 16384
    g = π * (3 - sqrt(5))
    [(sqrt(1 - z^2) * cos(g * k), sqrt(1 - z^2) * sin(g * k), z) for (k, z) in
     ((k, 1 - (2k + 1) / M) for k in 0:M-1)]
end
const _RAYS2 = [(cos(2π * (k + 0.5) / 4096), sin(2π * (k + 0.5) / 4096)) for k in 0:4095]

# planes n·x ≤ b of a convex body (cell-local), outward by its centroid; and the centroid
function _body_planes(P::_Poly3)
    sc = maximum(p -> maximum(abs, p), _vertices(P))
    v, f = _indexed(P, 1e-12 * sc)
    c = _polyhedron_centroid(v, f)
    planes = Tuple{NTuple{3,Float64},Float64}[]
    for face in f
        n = _face_normal(v, face)
        nn = _norm(n)
        nn > 0 || continue
        n = n ./ nn
        b = _dot(n, v[face[1]])
        _dot(n, c) > b && (n = .-n; b = -b)
        push!(planes, (n, b))
    end
    return planes, c
end
function _body_planes(P::_Poly2)
    c = _polygon_centroid(P.pts)
    planes = Tuple{NTuple{2,Float64},Float64}[]
    m = length(P.pts)
    for k in 1:m
        p, q = P.pts[k], P.pts[mod1(k + 1, m)]
        n = (q[2] - p[2], p[1] - q[1])
        nn = _norm(n)
        nn > 0 || continue
        n = n ./ nn
        b = _dot(n, p)
        _dot(n, c) > b && (n = .-n; b = -b)
        push!(planes, (n, b))
    end
    return planes, c
end

# ball of a pore (cell-local centre, radius); rounding in (0, 1]. The radius is moved by up to
# a few tenths of a percent away from tangency with any vertex, edge or face of the pore or of
# the cell (and their intersections): OCC splits nearly tangent curves differently on periodic
# copies, which Gmsh cannot match.
function _pore_ball(P::_Poly3, rounding, L)
    planes, c = _body_planes(P)
    R = (3 * _poly_volume(P) / 4π)^(1 / 3) / rounding
    return c, _clear_radius(R, _critical_distances(P, planes, c, L))
end
function _pore_ball(P::_Poly2, rounding, L)
    planes, c = _body_planes(P)
    R = sqrt(_poly_volume(P) / π) / rounding
    return c, _clear_radius(R, _critical_distances(P, planes, c, L))
end

_segdist(p, a, b) = (ab = b .- a; s = clamp(_dot(p .- a, ab) / _dot(ab, ab), 0.0, 1.0); _norm(p .- (a .+ s .* ab)))

function _edges(P::_Poly3)
    E = Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
    for f in P.faces, k in eachindex(f)
        push!(E, (f[k], f[mod1(k + 1, length(f))]))
    end
    return E
end
_edges(P::_Poly2) = [(P.pts[k], P.pts[mod1(k + 1, length(P.pts))]) for k in eachindex(P.pts)]

# distances from c at which a sphere centred at c is tangent to (or passes through) a feature
function _critical_distances(P, planes, c::NTuple{D,Float64}, L) where {D}
    d = Float64[]
    V = _bodyvertices(P)
    blo = ntuple(k -> minimum(v -> v[k], V), D)
    bhi = ntuple(k -> maximum(v -> v[k], V), D)
    tol = 1e-6 * maximum(L)
    # a point of a pore face (inside the pore's planes) and of the closed cell
    onpore(q) = all(((n, b),) -> _dot(n, q) <= b + tol, planes) && all(k -> blo[k] - tol <= q[k] <= bhi[k] + tol, 1:D)
    incell(q) = all(k -> -tol <= q[k] <= L[k] + tol, 1:D)
    for v in _bodyvertices(P)
        push!(d, _norm(v .- c))
    end
    E = _edges(P)
    for (a, b) in E
        push!(d, _segdist(c, a, b))
    end
    for (n, b) in planes
        push!(d, abs(b - _dot(n, c)))
    end
    for k in 1:D, x in (0.0, L[k])
        push!(d, abs(c[k] - x))                                        # cell face
        for (a, b) in E                                                # pore edge through the face
            (a[k] - x) * (b[k] - x) < 0 || continue
            q = a .+ ((x - a[k]) / (b[k] - a[k])) .* (b .- a)
            push!(d, _norm(q .- c))
        end
        for (n, b) in planes                                           # pore face ∩ cell face
            D == 3 || continue
            m = ntuple(j -> j == k ? 1.0 : 0.0, 3)
            u = _cross(n, m)
            nu = _norm(u)
            nu > 1e-9 || continue
            # point on both planes closest to c: c + α n + β m, n·x = b, x_k = x
            nm = n[k]
            rhs = [b - _dot(n, c), x - c[k]]
            abs(1 - nm^2) > 1e-12 || continue
            det = 1 - nm^2
            α = (rhs[1] - nm * rhs[2]) / det
            β = (rhs[2] - nm * rhs[1]) / det
            q = c .+ α .* n .+ β .* m
            onpore(q) && incell(q) && push!(d, _norm(q .- c))
        end
        for j in k+1:D, y in (0.0, L[j])                              # cell edge
            push!(d, sqrt((c[k] - x)^2 + (c[j] - y)^2))
            D == 3 || continue
            # pore face through the cell edge: the point on the edge line
            i = 6 - k - j
            for (nv, b) in planes
                abs(nv[i]) > 1e-9 || continue
                q = ntuple(l -> l == k ? x : l == j ? y : 0.0, 3)
                qi = (b - nv[k] * x - nv[j] * y) / nv[i]
                q = ntuple(l -> l == i ? qi : q[l], 3)
                onpore(q) && incell(q) && push!(d, _norm(q .- c))
            end
        end
    end
    for corner in Iterators.product(ntuple(k -> (0.0, L[k]), D)...)   # cell corners
        push!(d, _norm(corner .- c))
    end
    # beyond the farthest pore vertex the ball contains the pore: nothing to clear
    far = maximum(v -> _norm(v .- c), V)
    return filter(x -> x <= 1.01 * far, d)
end

# the radius nearest to R that keeps a relative margin δ from every critical distance, moving
# it by at most 0.5 %; if there is none, R (the generate retries handle the rare failures)
function _clear_radius(R, d; δ=0.002, maxmove=0.005)
    ok(r) = all(x -> abs(r - x) >= δ * x, d)
    ok(R) && return R
    best = R
    gap = maxmove * R
    for x in d, r in (x * (1 + δ) * (1 + 1e-9), x * (1 - δ) * (1 - 1e-9))
        abs(r - R) <= gap && ok(r) && ((best, gap) = (r, abs(r - R)))
    end
    return best
end

_bodyvertices(P::_Poly3) = _vertices(P)
_bodyvertices(P::_Poly2) = P.pts

# volume (area) of a convex body ∩ ball(c, R): exact when the ball contains the body
function _ball_volume(P, c::NTuple{D,Float64}, R) where {D}
    maximum(p -> _norm(p .- c), _bodyvertices(P)) <= R && return _poly_volume(P)
    planes, _ = _body_planes(P)
    rays = D == 3 ? _RAYS3 : _RAYS2
    s = 0.0
    @inbounds for u in rays
        t0, t1 = 0.0, R
        for (n, b) in planes
            a = _dot(n, u)
            h = b - _dot(n, c)
            if a > 0
                t1 = min(t1, h / a)
            elseif a < 0
                t0 = max(t0, h / a)
            elseif h < 0
                t1 = -1.0
            end
        end
        t1 > t0 && (s += t1^D - t0^D)
    end
    return D == 3 ? s * 4π / (3 * length(rays)) : s * π / length(rays)
end

# pore volume inside the cell of grain i with wall w and rounding r (0: polyhedral)
function _pore_volume(t::Tessellation{D}, i, w, r) where {D}
    P = _shrunk(t, i, w / 2)
    P === nothing && return 0.0
    r == 0 && return sum(_poly_volume, _cell_pieces(t, P); init=0.0)
    c, R = _pore_ball(P, r, t.cell.size)
    return sum(_ball_volume(Q, c .+ s .* t.cell.size, R) for (Q, s) in _cell_pieces(t, P; shifts=true); init=0.0)
end

_foam_density(t::Tessellation, w, r) = 1 - sum(_pore_volume(t, i, w, r) for i in eachindex(t.loops)) / prod(t.cell.size)

function _check_rounding(r, what)
    r isa Real && !(r isa Bool) && 0 <= r <= 1 || throw(ArgumentError("$what: rounding must be a number in [0, 1], got $(repr(r))"))
    return Float64(r)
end

# per grain, for a fixed wall: the pieces of the shrunk grain (with their planes and shifts),
# the centre and equivalent radius of its ball and the critical distances of the radius
function _pore_data(t::Tessellation{D}, w) where {D}
    L = t.cell.size
    out = []
    for i in eachindex(t.loops)
        P = _shrunk(t, i, w / 2)
        P === nothing && continue
        planes, c = _body_planes(P)
        r0 = D == 3 ? (3 * _poly_volume(P) / 4π)^(1 / 3) : sqrt(_poly_volume(P) / π)
        crit = _critical_distances(P, planes, c, L)
        pieces = [(Q, c .+ s .* L, first(_body_planes(Q)), _poly_volume(Q), maximum(p -> _norm(p .- (c .+ s .* L)), _bodyvertices(Q)))
                  for (Q, s) in _cell_pieces(t, P; shifts=true)]
        push!(out, (r0, crit, pieces))
    end
    return out
end

function _ball_volume(planes, c::NTuple{D,Float64}, R, vol, rmax) where {D}
    rmax <= R && return vol
    rays = D == 3 ? _RAYS3 : _RAYS2
    s = 0.0
    @inbounds for u in rays
        t0, t1 = 0.0, R
        for (n, b) in planes
            a = _dot(n, u)
            h = b - _dot(n, c)
            if a > 0
                t1 = min(t1, h / a)
            elseif a < 0
                t0 = max(t0, h / a)
            elseif h < 0
                t1 = -1.0
            end
        end
        t1 > t0 && (s += t1^D - t0^D)
    end
    return D == 3 ? s * 4π / (3 * length(rays)) : s * π / length(rays)
end

function _cached_density(data, r, V)
    pv = 0.0
    for (r0, crit, pieces) in data
        R = _clear_radius(r0 / r, crit)
        for (planes, cs, vol, rmax) in ((p[3], p[2], p[4], p[5]) for p in pieces)
            pv += _ball_volume(planes, cs, R, vol, rmax)
        end
    end
    return 1 - pv / V
end

"""
    pore_rounding(t, density; wall) -> rounding

Rounding of the closed-cell foam [`pores`](@ref)`(t; wall, rounding)` with the given relative
density (volume fraction of the walls) and wall thickness: the pores get rounder (thicker
nodes and edges) until the solid reaches `density`. Throws an error when the density is out of
reach: below the polyhedral foam with that wall (`rounding = 0`, reduce `wall`) or above
`rounding = 1` (increase `wall`).

The density is reached to a few 1e-3 (relative) with a few tens of pores, up to about 2e-2
with only a few pores: each ball radius is moved by up to 0.5 % away from tangencies that Gmsh cannot mesh,
so the density jumps slightly with the rounding. `volume_fraction` of the meshed RVE gives the
exact value.
"""
function pore_rounding(t::Tessellation{D}, density::Real; wall) where {D}
    0 < density < 1 || throw(ArgumentError("pore_rounding: density must be in (0, 1)"))
    w = Float64(wall)
    w >= 0 || throw(ArgumentError("pore_rounding: wall must be >= 0"))
    V = prod(t.cell.size)
    data = _pore_data(t, w)
    ρ(r) = _cached_density(data, r, V)
    ρ0, ρ1 = _foam_density(t, w, 0.0), ρ(1.0)
    density < ρ0 && throw(ArgumentError(@sprintf("pore_rounding: the density %.4g is below the polyhedral foam with wall = %.4g (%.4g); reduce wall",
                                                 density, w, ρ0)))
    density > ρ1 && throw(ArgumentError(@sprintf("pore_rounding: the density %.4g is above the fully rounded foam (rounding = 1) with wall = %.4g (%.4g); increase wall",
                                                 density, w, ρ1)))
    lo, hi = 1e-3, 1.0
    ρlo, ρhi = ρ(lo), ρ1
    for _ in 1:40
        m = (lo + hi) / 2
        ρm = ρ(m)
        ρm < density ? ((lo, ρlo) = (m, ρm)) : ((hi, ρhi) = (m, ρm))
        hi - lo < 1e-8 && break
    end
    # the density may jump inside the last bracket: the end closer to the target
    return abs(ρlo - density) <= abs(ρhi - density) ? lo : hi
end

# OCC parametrises a sphere with poles along its z axis and the seam (u = 0) towards +x; a
# plane (cell face or pore face) that passes near a pole or touches the seam tangentially
# leaves curves split differently on periodic copies, or overlapping facets. Among a fixed
# set of rotations, take the one whose poles and seam stay farthest from every plane that
# cuts the ball.
const _BALL_ROTATIONS = let
    out = Tuple{NTuple{3,Float64},Float64}[((0.0, 0.0, 1.0), 0.0)]
    g = π * (3 - sqrt(5))
    for k in 0:22
        z = 1 - (2k + 1) / 23
        push!(out, ((sqrt(1 - z^2) * cos(g * k), sqrt(1 - z^2) * sin(g * k), z), 0.4 + 2.3 * mod(k * 0.618034, 1)))
    end
    out
end

function _rodrigues(axis, θ, v)
    c, s = cos(θ), sin(θ)
    return v .* c .+ _cross(axis, v) .* s .+ axis .* (_dot(axis, v) * (1 - c))
end

function _ball_rotation(cs::NTuple{3,Float64}, R, planes)
    cuts = [(n, b) for (n, b) in planes if abs(b - _dot(n, cs)) < R]
    isempty(cuts) && return first(_BALL_ROTATIONS)
    best, bestm = first(_BALL_ROTATIONS), -Inf
    for (axis, θ) in _BALL_ROTATIONS
        a = _rodrigues(axis, θ, (1.0, 0.0, 0.0))
        e = _rodrigues(axis, θ, (0.0, 0.0, 1.0))
        m = Inf
        for (n, b) in cuts
            h = (_dot(n, cs) - b) / R
            na, ne = _dot(n, a), _dot(n, e)
            m = min(m, abs(h + ne), abs(h - ne))                      # poles
            for k in 0:64                                             # seam: no near-tangency
                φ = π * k / 64
                g = h + sin(φ) * na + cos(φ) * ne
                dg = cos(φ) * na - sin(φ) * ne
                m = min(m, max(abs(g), abs(dg)))
            end
        end
        m > bestm && ((best, bestm) = ((axis, θ), m))
    end
    return best
end

function _disk_rotation(cs::NTuple{2,Float64}, R, planes)
    cuts = [(n, b) for (n, b) in planes if abs(b - _dot(n, cs)) < R]
    isempty(cuts) && return 0.0
    best, bestm = 0.0, -Inf
    for k in 0:23
        θ = 2π * k / 24 + (k > 0 ? 0.1 : 0.0)
        p = cs .+ R .* (cos(θ), sin(θ))
        m = minimum(abs(_dot(n, p) - b) for (n, b) in cuts) / R
        m > bestm && ((best, bestm) = (θ, m))
    end
    return best
end
