# Random short fibres: non-overlapping capsules of finite length with prescribed orientations,
# placed by random sequential addition (orientation fixed per fibre, only the position redrawn).

const _V3 = NTuple{3,Float64}
_dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

# squared distance between the segments p1-q1 and p2-q2 (Ericson, Real-Time Collision Detection, 5.1.9)
function _segdist2(p1::_V3, q1::_V3, p2::_V3, q2::_V3)
    d1 = q1 .- p1; d2 = q2 .- p2; r = p1 .- p2
    a = _dot3(d1, d1); e = _dot3(d2, d2); f = _dot3(d2, r)
    s = 0.0; t = 0.0
    if a <= 0 && e <= 0
        return _dot3(r, r)
    elseif a <= 0
        t = clamp(f / e, 0.0, 1.0)
    else
        c = _dot3(d1, r)
        if e <= 0
            s = clamp(-c / a, 0.0, 1.0)
        else
            b = _dot3(d1, d2); den = a * e - b * b
            s = den > 1e-14 * a * e ? clamp((b * f - c * e) / den, 0.0, 1.0) : 0.0
            t = (b * s + f) / e
            if t < 0
                t = 0.0; s = clamp(-c / a, 0.0, 1.0)
            elseif t > 1
                t = 1.0; s = clamp((b - c) / a, 0.0, 1.0)
            end
        end
    end
    w = (p1 .+ d1 .* s) .- (p2 .+ d2 .* t)
    return _dot3(w, w)
end

# segments with a radius and all their periodic images that reach the cell, binned on a grid
struct _SegGrid
    L::_V3
    per::NTuple{3,Bool}
    nb::NTuple{3,Int}
    reach::Float64
    bins::Dict{NTuple{3,Int},Vector{Int}}
    A::Vector{_V3}
    B::Vector{_V3}
    R::Vector{Float64}
end
_SegGrid(L, per, reach, w) = _SegGrid(L, per, ntuple(d -> clamp(floor(Int, L[d] / w), 1, 64), 3), reach,
                                      Dict{NTuple{3,Int},Vector{Int}}(), _V3[], _V3[], Float64[])
_bin(g::_SegGrid, x, d) = clamp(floor(Int, x / g.L[d] * g.nb[d]), 0, g.nb[d] - 1)

function _push_seg!(g::_SegGrid, a::_V3, b::_V3, r; skip=0)
    lo = min.(a, b) .- (r + g.reach); hi = max.(a, b) .+ (r + g.reach)
    ks = ntuple(d -> g.per[d] && d != skip ? (floor(Int, -hi[d] / g.L[d]):ceil(Int, (g.L[d] - lo[d]) / g.L[d])) : (0:0), 3)
    for k in Iterators.product(ks...)
        s = k .* g.L
        l2 = lo .+ s; h2 = hi .+ s
        all(d -> h2[d] > 0 && l2[d] < g.L[d], 1:3) || continue
        push!(g.A, a .+ s); push!(g.B, b .+ s); push!(g.R, r)
        id = length(g.R)
        for c in Iterators.product(ntuple(d -> _bin(g, l2[d], d):_bin(g, h2[d], d), 3)...)
            push!(get!(g.bins, c, Int[]), id)
        end
    end
end

function _seg_clear(g::_SegGrid, a::_V3, b::_V3, r, gap)
    lo = min.(a, b); hi = max.(a, b)
    for c in Iterators.product(ntuple(d -> _bin(g, lo[d], d):_bin(g, hi[d], d), 3)...)
        for j in get(g.bins, c, ())
            _segdist2(a, b, g.A[j], g.B[j]) < (r + g.R[j] + gap)^2 && return false
        end
    end
    return true
end

# the candidate and its periodic copies that reach the cell keep clear of everything stored
function _seg_clear_all(g::_SegGrid, a, b, r, gap)
    lo = min.(a, b) .- r; hi = max.(a, b) .+ r
    ks = ntuple(d -> g.per[d] ? (floor(Int, -hi[d] / g.L[d]):ceil(Int, (g.L[d] - lo[d]) / g.L[d])) : (0:0), 3)
    for k in Iterators.product(ks...)
        s = k .* g.L
        all(d -> hi[d] + s[d] > 0 && lo[d] + s[d] < g.L[d], 1:3) || continue
        _seg_clear(g, a .+ s, b .+ s, r, gap) || return false
    end
    return true
end

# a capsule longer than the cell must also keep clear of its own images
function _self_ok(L, per, a, b, r, gap)
    ext = abs.(b .- a) .+ (2r + gap)
    any(d -> per[d] && ext[d] > L[d], 1:3) || return true
    K = ntuple(d -> per[d] ? ceil(Int, ext[d] / L[d]) : 0, 3)
    for k in Iterators.product(ntuple(d -> -K[d]:K[d], 3)...)
        all(==(0), k) && continue
        s = k .* L
        _segdist2(a, b, a .+ s, b .+ s) < (2r + gap)^2 && return false
    end
    return true
end

# bands to the faces, edges and corners (slivers and Gmsh periodicity failures): the surface
# nearly tangent to a face plane, the cylinder-cap junction circles and the poles of the cap
# spheres (as build!(::Capsule) makes them) near a face, an edge or corner nearly tangent
function _capsule_bands_ok(L, per, a::_V3, b::_V3, r, f)
    f == 0 && return true
    u = (b .- a) ./ sqrt(_dot3(b .- a, b .- a))
    for d in 1:3
        e = r * sqrt(max(0.0, 1 - u[d]^2))
        vals = (min(a[d], b[d]) - r, max(a[d], b[d]) + r, a[d] - e, a[d] + e, b[d] - e, b[d] + e)
        # cap poles along z (build! makes the spheres there) and the seam / pole heights at the centres
        vals = d == 3 ? (vals..., a[d] - r, a[d] + r, b[d] - r, b[d] + r, a[d], b[d]) : (vals..., a[d], b[d])
        for v in vals
            m = per[d] ? abs(v - L[d] * round(v / L[d])) : min(abs(v), abs(v - L[d]))
            m < f && return false
        end
    end
    # cell edges and corners (all lattice copies that the capsule can reach)
    kr = ntuple(d -> floor(Int, (min(a[d], b[d]) - r - f) / L[d]):ceil(Int, (max(a[d], b[d]) + r + f) / L[d]), 3)
    for d in 1:3
        t = d == 1 ? (2, 3) : d == 2 ? (1, 3) : (1, 2)
        for ka in kr[t[1]], kb in kr[t[2]]
            (!per[t[1]] && !(ka in (0, 1))) && continue
            (!per[t[2]] && !(kb in (0, 1))) && continue
            p = ntuple(i -> i == d ? (kr[d][1] - 1) * L[d] : i == t[1] ? ka * L[i] : kb * L[i], 3)
            q = ntuple(i -> i == d ? (kr[d][end] + 1) * L[d] : p[i], 3)
            abs(sqrt(_segdist2(a, b, p, q)) - r) < f && return false
        end
    end
    for k in Iterators.product(kr...)
        any(d -> !per[d] && !(k[d] in (0, 1)), 1:3) && continue
        p = Float64.(k .* L)
        abs(sqrt(_segdist2(a, b, p, p)) - r) < f && return false
    end
    return true
end

_capvol(l, r) = π * r^2 * l + 4π * r^3 / 3
_iso_dir(rng) = (z = 2rand(rng) - 1; φ = 2π * rand(rng); s = sqrt(1 - z^2); (s * cos(φ), s * sin(φ), z))

function _orientation(orientation, axis, what)
    axis isa Integer && !(axis isa Bool) && axis in 1:3 || throw(ArgumentError("$what: axis must be 1, 2 or 3, got $(repr(axis))"))
    if orientation === :isotropic
        return _iso_dir
    elseif orientation === :planar
        t = Tuple(d for d in 1:3 if d != axis)
        return rng -> (φ = 2π * rand(rng); ntuple(d -> d == t[1] ? cos(φ) : d == t[2] ? sin(φ) : 0.0, 3))
    elseif orientation isa Union{Tuple,AbstractVector} && length(orientation) == 3 && all(x -> x isa Real, orientation)
        u = _normalize(_vec(Val(3), orientation, "$what orientation"), "$what orientation")
        return _ -> u
    elseif orientation isa Symbol || orientation isa AbstractString || orientation isa Union{Tuple,AbstractVector,AbstractMatrix}
        throw(ArgumentError("$what: orientation must be :isotropic, :planar, a direction (3 numbers) or a function rng -> direction, got $(repr(orientation))"))
    end
    return rng -> _normalize(_vec(Val(3), orientation(rng), "$what sampled direction"), "$what sampled direction")
end

function _lengths(rng, V, len, r, fraction, n, what)
    fraction === nothing || (fraction isa Real && 0 < fraction < 1) ||
        throw(ArgumentError("$what: fraction must be a number in (0, 1), got $(repr(fraction))"))
    n === nothing || (n isa Integer && !(n isa Bool) && n >= 1) ||
        throw(ArgumentError("$what: n must be a positive integer, got $(repr(n))"))
    if len isa AbstractVector
        (fraction === nothing && n === nothing) || throw(ArgumentError("$what: with a vector of lengths, give neither fraction nor n"))
        isempty(len) && throw(ArgumentError("$what: empty vector of lengths"))
        l = [_positive(x, "$what length") for x in len]
    elseif len isa Real
        (fraction === nothing) != (n === nothing) || throw(ArgumentError("$what: give exactly one of fraction and n"))
        l0 = _positive(len, "$what length")
        l = fill(l0, n === nothing ? max(1, round(Int, fraction * V / _capvol(l0, r))) : Int(n))
    elseif len === nothing || len isa Union{AbstractString,Symbol}
        throw(ArgumentError("$what: length must be a number, a vector of lengths or a function rng -> length"))
    else
        (fraction === nothing) != (n === nothing) || throw(ArgumentError("$what: with a length sampler give fraction or n"))
        l = Float64[]
        if n !== nothing
            for _ in 1:n
                push!(l, _positive(len(rng), "$what sampled length"))
            end
        else
            vol = 0.0
            while true
                x = _positive(len(rng), "$what sampled length")
                vol + _capvol(x, r) > fraction * V && break
                push!(l, x); vol += _capvol(x, r)
                length(l) > 10^6 && throw(ArgumentError("$what: more than 10^6 fibres; check the sampler"))
            end
            isempty(l) && throw(ArgumentError("$what: the first sampled fibre is larger than the requested fraction"))
        end
    end
    return sort!(l; rev=true)
end

"""
    random_capsules(cell; length, radius, fraction, n, orientation = :isotropic, axis = 3,
                    min_gap = :auto, face_gap = :auto, boundary = :periodic, avoid = nothing,
                    rng = Random.default_rng(), max_tries = 50_000, strict = true) -> Vector{Capsule}

Random short fibres: non-overlapping capsules (cylinders with hemispherical ends) of radius
`radius` placed by random sequential addition. `length` is the distance between the centres of
the two caps (total length `length + 2 radius`): a number with `fraction` or `n` (exactly one),
a vector of lengths, or a sampler `rng -> length` with `fraction` or `n`.

- `orientation`: `:isotropic` (uniform on the sphere), `:planar` (uniform in the plane normal to
  `axis`), a direction `(x, y, z)` (all aligned with it) or a function `rng -> direction`. Each
  fibre draws its orientation once; failed attempts only move it, so the orientations are not
  biased by the rejections. Orientations in which the fibre cannot exist at all (longer than
  the cell and touching its own periodic images, or not fitting with `boundary = :inside`) are
  redrawn, which biases the distribution only for fibres longer than the cell.
  [`orientation_tensor`](@ref) gives their second moment.
- The fibres are returned from the longest to the shortest.
- `min_gap` (default `0.1 radius`), `face_gap` (default `min_gap`), `boundary`, `avoid`, `rng`,
  `max_tries` and `strict` as in [`random_spheres`](@ref). The bands to the faces reject fibres
  whose surface, cap junctions or cap poles nearly touch a face, edge or corner.
- `fraction` is the volume fraction of whole capsules (fibres crossing free faces are cut).

Random sequential addition of equal capsules with the default gap reaches about 0.26, 0.19,
0.15 and 0.09 (isotropic) and 0.33 to 0.31 (aligned) for total aspect ratios
`(length + 2r) / 2r` of 2, 5, 10 and 20; above that an error reports what was placed. Put the
result in a phase without `fuse` and refine its interface to `radius / 2`–`radius` (a coarse
mesh of thin fibres can make Gmsh fail):
`Phase(:fibres, random_capsules(cell; length = 0.3, radius = 0.04, n = 10); refine = 0.02)`.
"""
function random_capsules(cell::Cell{3}; length=nothing, radius=nothing, fraction=nothing, n=nothing, orientation=:isotropic,
                         axis=3, min_gap=:auto, face_gap=:auto, boundary=:periodic, avoid=nothing,
                         rng=Random.default_rng(), max_tries=50_000, strict=true)
    what = "random_capsules"
    rng = _rng(rng)
    r = _positive(radius === nothing ? throw(ArgumentError("$what: give the radius")) : radius, "$what radius")
    boundary in (:periodic, :inside) || throw(ArgumentError("$what: boundary must be :periodic or :inside, got $(repr(boundary))"))
    max_tries isa Integer && max_tries >= 1 || throw(ArgumentError("$what: max_tries must be a positive integer, got $(repr(max_tries))"))
    gap = min_gap === :auto ? 0.1r :
          min_gap isa Real && isfinite(min_gap) && min_gap >= 0 ? Float64(min_gap) :
          throw(ArgumentError("$what: min_gap must be :auto or a finite number >= 0, got $(repr(min_gap))"))
    fgap = face_gap === :auto ? gap :
           face_gap isa Real && isfinite(face_gap) && face_gap >= 0 ? Float64(face_gap) :
           throw(ArgumentError("$what: face_gap must be :auto or a finite number >= 0, got $(repr(face_gap))"))
    orient = _orientation(orientation, axis, what)
    L = Tuple(Float64.(cell.size))
    per = cell.periodic
    inside = boundary === :inside
    lens = _lengths(rng, volume(cell), length, r, fraction, n, what)
    all(d -> 2r + (inside ? 2fgap : per[d] ? gap : 0.0) < L[d], 1:3) ||
        throw(ArgumentError("$what: the fibre diameter does not fit in the cell"))
    # an orientation is impossible when the fibre cannot fit inside (boundary = :inside) or would
    # touch its own periodic images wherever it is placed (fibres longer than the cell)
    fits(u, l) = (!inside || all(d -> abs(u[d]) * l + 2r + 2fgap < L[d], 1:3)) &&
                 _self_ok(L, per, (0.0, 0.0, 0.0), u .* l, r, gap)
    if orientation isa Union{Tuple,AbstractVector}
        u = orient(rng)
        fits(u, lens[1]) ||
            throw(ArgumentError("$what: the longest fibre in this direction does not fit in the cell " *
                                (inside ? "(boundary = :inside)" : "(it touches its own periodic images)") *
                                (any(d -> per[d] && abs(u[d]) > 1 - 1e-12, 1:3) ? "; for fibres spanning the cell use random_fibres" : "")))
    end
    grid = _SegGrid(L, per, r + gap, 2r + gap)
    o = Tuple(Float64.(cell.origin))
    for g in _flat_items(avoid)
        g isa Geometry{3} || throw(ArgumentError("avoid: a 3D packing cannot avoid the $(nameof(typeof(g)))"))
        a, b, ro, fax = _segment(g, (0.0, 0.0, 0.0))
        a, b = Tuple(Float64.(a .- o)), Tuple(Float64.(b .- o))
        if fax > 0                               # infinite fibre: a segment longer than any candidate reach
            m = lens[1] + 2r + gap + ro
            a = ntuple(d -> d == fax ? -m : a[d], 3); b = ntuple(d -> d == fax ? L[d] + m : b[d], 3)
        end
        _push_seg!(grid, a, b, ro; skip=fax)
    end
    out = Capsule[]
    for l in lens
        u = orient(rng)
        redraws = 0
        while !fits(u, l)
            (redraws += 1) > 1000 &&
                throw(ArgumentError("$what: no sampled orientation lets a fibre of length $l fit in the cell; shorten the fibres"))
            u = orient(rng)
        end
        placed = false
        for _ in 1:max_tries
            c = ntuple(d -> rand(rng) * L[d], 3)
            a = c .- u .* (l / 2); b = c .+ u .* (l / 2)
            if inside
                all(d -> min(a[d], b[d]) - r >= fgap && max(a[d], b[d]) + r <= L[d] - fgap, 1:3) || continue
            else
                _capsule_bands_ok(L, per, a, b, r, fgap) || continue
            end
            _self_ok(L, per, a, b, r, gap) || continue
            _seg_clear_all(grid, a, b, r, gap) || continue
            _push_seg!(grid, a, b, r)
            push!(out, Capsule(o .+ a, o .+ b, r))
            placed = true
            break
        end
        placed || break
    end
    if Base.length(out) < Base.length(lens)
        V = volume(cell)
        got = sum(_capvol(l, r) for l in lens[1:Base.length(out)]; init=0.0) / V
        want = sum(_capvol(l, r) for l in lens) / V
        msg = @sprintf("%s: placed %d of %d fibres (fraction %.3f of the requested %.3f) after %d attempts for the next one; random sequential addition of capsules saturates early, more so for slender isotropic fibres (see the docstring); reduce fraction, length, min_gap or face_gap",
                       what, Base.length(out), Base.length(lens), got, want, max_tries)
        avoid === nothing || (msg *= "; the inclusions to avoid take part of the cell")
        strict ? error(msg) : @warn(msg)
    end
    return out
end

"""
    orientation_tensor(items) -> Matrix{Float64}

Second-order orientation tensor `A = (1/N) Σ pᵢ pᵢᵀ` (Advani-Tucker) of capsules, cylinders or
fibres (or translates of these), with `pᵢ` the unit axis of each item. `A ≈ I/3` for isotropic
fibres, `diag(1/2, 1/2, 0)` for random orientations in the xy plane.
"""
function orientation_tensor(items)
    A = zeros(3, 3)
    N = 0
    for it0 in items
        it = it0
        while it isa Translate
            it = it.child
        end
        v = it isa Capsule ? it.p2 .- it.p1 :
            it isa Cylinder ? it.axis :
            it isa Fibre ? it.direction :
            throw(ArgumentError("orientation_tensor: items must be capsules, cylinders or fibres, got $(nameof(typeof(it)))"))
        l = _norm(v)
        l > 0 || continue
        p = v ./ l
        for i in 1:3, j in 1:3
            A[i, j] += p[i] * p[j]
        end
        N += 1
    end
    N > 0 || throw(ArgumentError("orientation_tensor: no items"))
    return A ./ N
end
