# Random packings of non-overlapping spheres, disks and aligned fibres (random sequential
# addition with gaps to the faces, edges and corners of the cell).

const _SAT = Dict(1 => 0.7476, 2 => 0.547, 3 => 0.384)   # RSA saturation of equal segments / disks / spheres

_ballvol(::Val{3}, r) = 4π * r^3 / 3
_ballvol(::Val{2}, r) = π * r^2

_rng(rng::AbstractRNG) = rng
_rng(seed::Integer) = seed >= 0 ? Random.Xoshiro(seed) : throw(ArgumentError("rng: a seed must be >= 0, got $seed"))
_rng(x) = throw(ArgumentError("rng must be an AbstractRNG or an integer seed, got $(typeof(x))"))

"""
Radii (decreasing) of the packing from the keyword combination; `V` is the volume (area in 2D)
of the region, `what` names the caller in errors.
"""
function _radii(::Val{D}, rng, V, radius, fraction, n, what) where {D}
    fraction === nothing || fraction isa Real ||
        throw(ArgumentError("$what: fraction must be a number in (0, 1), got $(repr(fraction))"))
    n === nothing || (n isa Integer && !(n isa Bool)) ||
        throw(ArgumentError("$what: n must be an integer, got $(repr(n))"))
    fraction === nothing || 0 < fraction < 1 ||
        throw(ArgumentError("$what: fraction must be in (0, 1), got $fraction"))
    n === nothing || n >= 1 || throw(ArgumentError("$what: n must be at least 1, got $n"))
    if radius isa AbstractVector
        (fraction === nothing && n === nothing) ||
            throw(ArgumentError("$what: with a vector of radii, give neither fraction nor n"))
        isempty(radius) && throw(ArgumentError("$what: empty vector of radii"))
        r = [_positive(x, "$what radius") for x in radius]
    elseif radius isa Real
        (fraction === nothing) != (n === nothing) ||
            throw(ArgumentError("$what: give exactly two of radius, fraction and n"))
        r0 = _positive(radius, "$what radius")
        m = n === nothing ? max(1, round(Int, fraction * V / _ballvol(Val(D), r0))) : Int(n)
        r = fill(r0, m)
    elseif radius === nothing
        (fraction !== nothing && n !== nothing) ||
            throw(ArgumentError("$what: give exactly two of radius, fraction and n"))
        r = fill((fraction * V / (n * _ballvol(Val(D), 1.0)))^(1 / D), Int(n))
    elseif !(radius isa AbstractString || radius isa Symbol)          # a sampler: any callable rng -> r
        fraction !== nothing || n !== nothing ||
            throw(ArgumentError("$what: a radius sampler needs fraction or n"))
        (fraction === nothing) != (n === nothing) ||
            throw(ArgumentError("$what: with a radius sampler give fraction or n, not both"))
        r = Float64[]
        if n !== nothing
            for _ in 1:n
                push!(r, _positive(radius(rng), "$what sampled radius"))
            end
        else
            vol = 0.0
            while true
                x = _positive(radius(rng), "$what sampled radius")
                vol + _ballvol(Val(D), x) > fraction * V && break
                push!(r, x); vol += _ballvol(Val(D), x)
                length(r) > 10^6 && throw(ArgumentError("$what: more than 10^6 inclusions; check the sampler"))
            end
            isempty(r) && throw(ArgumentError("$what: the first sampled inclusion is larger than the requested fraction"))
        end
    else
        throw(ArgumentError("$what: radius must be a number, a vector of radii or a function rng -> radius"))
    end
    return sort!(r; rev=true)
end

"""
Random sequential addition in the box `[0, L]` (cell-local coordinates). Returns the centres
(in the order of `radii`) and the number placed.
"""
function _rsa(rng, L::NTuple{D,Float64}, per::NTuple{D,Bool}, radii, gap, fgap, inside::Bool, maxtries;
              obst=nothing, propose=nothing) where {D}
    rmax = radii[1]
    w = 2rmax + gap
    nb = ntuple(d -> max(1, floor(Int, L[d] / w)), D)
    bins = Dict{NTuple{D,Int},Vector{Int}}()
    binof(x) = ntuple(d -> clamp(floor(Int, x[d] / L[d] * nb[d]), 0, nb[d] - 1), D)
    X = NTuple{D,Float64}[]
    R = Float64[]
    function neighbours(b)
        rng_d = ntuple(d -> nb[d] < 3 ? (0:nb[d]-1) : (b[d]-1:b[d]+1), D)
        out = NTuple{D,Int}[]
        for k in Iterators.product(rng_d...)
            kk = ntuple(D) do d
                per[d] ? mod(k[d], nb[d]) : k[d]
            end
            all(d -> 0 <= kk[d] < nb[d], 1:D) && push!(out, kk)
        end
        return unique!(out)
    end
    function face_ok(x, r)
        fgap == 0 && return true
        s = ntuple(d -> min(x[d], L[d] - x[d]), D)
        for d in 1:D
            (abs(s[d] - r) < fgap || s[d] < fgap) && return false
        end
        if D == 3
            for (a, b) in ((1, 2), (1, 3), (2, 3))
                abs(hypot(s[a], s[b]) - r) < fgap && return false
            end
        end
        abs(sqrt(sum(abs2, s)) - r) < fgap && return false
        return true
    end
    for r in radii
        lo = ntuple(d -> inside ? r + fgap : 0.0, D)
        hi = ntuple(d -> inside ? L[d] - r - fgap : L[d], D)
        all(d -> hi[d] > lo[d], 1:D) || break
        placed = false
        for _ in 1:maxtries
            x = propose === nothing ? ntuple(d -> lo[d] + rand(rng) * (hi[d] - lo[d]), D) : propose(rng, lo, hi)
            x === nothing && continue
            inside || face_ok(x, r) || continue
            ok = true
            for b in neighbours(binof(x)), j in get(bins, b, Int[])
                dist2 = 0.0
                for d in 1:D
                    δ = x[d] - X[j][d]
                    per[d] && (δ -= L[d] * round(δ / L[d]))
                    dist2 += δ^2
                end
                if dist2 < (r + R[j] + gap)^2
                    ok = false
                    break
                end
            end
            ok || continue
            obst === nothing || _clear(obst, x, r + gap) || continue
            push!(X, x); push!(R, r)
            push!(get!(bins, binof(x), Int[]), length(X))
            placed = true
            break
        end
        placed || break
    end
    return X
end


# ---------------------------------------------------------------- obstacles (avoid) and clusters

# obstacles as segments with a radius (spheres: a == b), in cell-local coordinates, with all the
# periodic images that a candidate inside the cell can reach
struct _Obstacles{D}
    a::Vector{NTuple{D,Float64}}
    b::Vector{NTuple{D,Float64}}
    r::Vector{Float64}
end

_flat_items(x::Geometry) = Geometry[x]
_flat_items(::Nothing) = Geometry[]
function _flat_items(x::Union{AbstractVector,Tuple})
    out = Geometry[]
    for el in x
        if el isa Geometry
            push!(out, el)
        elseif el isa Union{AbstractVector,Tuple}
            for g in el
                g isa Geometry || throw(ArgumentError("avoid: expected geometries, got $(typeof(g))"))
                push!(out, g)
            end
        else
            throw(ArgumentError("avoid: expected geometries or collections of geometries, got $(typeof(el))"))
        end
    end
    return out
end
_flat_items(x) = throw(ArgumentError("avoid: expected a geometry or a collection of geometries, got $(typeof(x))"))

const _AVOID_TYPES = "Sphere, Disk, Capsule, Cylinder, Ellipsoid, Ellipse (by the capsule that contains them), " *
                     "Fibre parallel to a cell axis, or a Translate of these"

# (a, b, radius, axis of an infinite fibre or 0) of a geometry, in global coordinates
_segment(g::Union{Sphere,Disk}, s) = (g.center .+ s, g.center .+ s, g.radius, 0)
_segment(g::Capsule, s) = (g.p1 .+ s, g.p2 .+ s, g.radius, 0)
_segment(g::Cylinder, s) = (g.base .+ s, g.base .+ g.axis .+ s, g.radius, 0)
function _segment(g::Union{Ellipsoid,Ellipse}, s)
    k = argmax(collect(g.radii))
    b = maximum(g.radii[i] for i in eachindex(g.radii) if i != k)
    e = ntuple(i -> i == k ? g.radii[k] - b : 0.0, length(g.radii))
    return (g.center .+ s .- e, g.center .+ s .+ e, b, 0)
end
function _segment(g::Fibre, s)
    k = findfirst(x -> abs(abs(x) - 1) < 1e-12, g.direction)
    k === nothing && throw(ArgumentError("avoid: only fibres parallel to a cell axis are supported"))
    e = ntuple(i -> i == k ? 1.0 : 0.0, 3)
    return (g.point .+ s .- 1e6 .* e, g.point .+ s .+ 1e6 .* e, g.radius, k)
end
_segment(g::Translate, s) = _segment(g.child, s .+ g.shift)
_segment(g, s) = throw(ArgumentError("avoid: unsupported geometry $(nameof(typeof(g))); supported: $_AVOID_TYPES"))

function _obstacles(::Val{D}, avoid, origin, L, per, reach) where {D}
    items = _flat_items(avoid)
    isempty(items) && return nothing
    A, B, Rs = NTuple{D,Float64}[], NTuple{D,Float64}[], Float64[]
    for g in items
        g isa Geometry{D} || throw(ArgumentError("avoid: a $(D)D packing cannot avoid the $(nameof(typeof(g))) " *
                                                 "of dimension $(3 - (D == 3))"))
        a, b, r, fax = _segment(g, ntuple(_ -> 0.0, D))
        a, b = Tuple(Float64.(a .- origin)), Tuple(Float64.(b .- origin))
        m = r + reach
        ks = ntuple(D) do d
            (per[d] && d != fax) || return 0:0
            lo, hi = min(a[d], b[d]) - m, max(a[d], b[d]) + m
            ceil(Int, -hi / L[d]):floor(Int, (L[d] - lo) / L[d])
        end
        for k in Iterators.product(ks...)
            sh = ntuple(d -> k[d] * L[d], D)
            push!(A, a .+ sh); push!(B, b .+ sh); push!(Rs, r)
        end
    end
    return _Obstacles{D}(A, B, Rs)
end

# true if the ball (x, ρ) keeps clear of every obstacle (ρ = radius + gap)
function _clear(o::_Obstacles{D}, x, ρ) where {D}
    @inbounds for i in eachindex(o.r)
        a, b, lim = o.a[i], o.b[i], o.r[i] + ρ
        any(d -> x[d] < min(a[d], b[d]) - lim || x[d] > max(a[d], b[d]) + lim, 1:D) && continue
        ab = b .- a
        l2 = sum(abs2, ab)
        t = l2 == 0 ? 0.0 : clamp(sum((x .- a) .* ab) / l2, 0.0, 1.0)
        sum(abs2, x .- (a .+ t .* ab)) < lim^2 && return false
    end
    return true
end

# cluster = (n | centers, radius, share = 1): proposal of candidate centres
function _cluster(::Val{D}, cluster, rng, origin, L, per, what) where {D}
    cluster isa NamedTuple || throw(ArgumentError("$what: cluster must be a named tuple (n = ..., radius = ...), got $(repr(cluster))"))
    bad = setdiff(keys(cluster), (:n, :radius, :share, :centers))
    isempty(bad) || throw(ArgumentError("$what: unknown cluster keys $(join(bad, ", ")); use n, radius, share, centers"))
    haskey(cluster, :n) == haskey(cluster, :centers) &&
        throw(ArgumentError("$what: cluster needs exactly one of n and centers"))
    R = get(cluster, :radius, nothing)
    R isa Real && R > 0 || throw(ArgumentError("$what: cluster radius must be a positive number"))
    share = get(cluster, :share, 1.0)
    share isa Real && 0 < share <= 1 || throw(ArgumentError("$what: cluster share must be in (0, 1], got $(repr(share))"))
    if haskey(cluster, :n)
        n = cluster.n
        n isa Integer && !(n isa Bool) && n >= 1 || throw(ArgumentError("$what: cluster n must be a positive integer, got $(repr(n))"))
        C = [ntuple(d -> rand(rng) * L[d], D) for _ in 1:n]
    else
        cs = cluster.centers
        cs isa Union{AbstractVector,Tuple} && !isempty(cs) || throw(ArgumentError("$what: cluster centers must be a non-empty list of points"))
        C = [Tuple(Float64.(_vec(Val(D), c, "cluster center") .- origin)) for c in cs]
    end
    function propose(rng, lo, hi)
        if share < 1 && rand(rng) >= share
            return ntuple(d -> lo[d] + rand(rng) * (hi[d] - lo[d]), D)
        end
        c = C[rand(rng, 1:length(C))]
        v = ntuple(_ -> randn(rng), D)
        nv = sqrt(sum(abs2, v))
        ρ = R * rand(rng)^(1 / D) / nv
        x = ntuple(d -> per[d] ? mod(c[d] + ρ * v[d], L[d]) : c[d] + ρ * v[d], D)
        all(d -> lo[d] <= x[d] <= hi[d], 1:D) || return nothing
        return x
    end
    return propose, length(C), Float64(R), Float64(share)
end

function _packing(::Val{D}, what, L, per, V, radius, fraction, n, min_gap, face_gap, boundary, rng, max_tries, strict;
                  avoid=nothing, cluster=nothing, origin=ntuple(_ -> 0.0, D)) where {D}
    rng = _rng(rng)
    boundary in (:periodic, :inside) || throw(ArgumentError("$what: boundary must be :periodic or :inside, got $(repr(boundary))"))
    max_tries isa Integer && max_tries >= 1 || throw(ArgumentError("$what: max_tries must be a positive integer, got $(repr(max_tries))"))
    radii = _radii(Val(D), rng, V, radius, fraction, n, what)
    gap = min_gap === :auto ? 0.1 * radii[end] :
          min_gap isa Real && isfinite(min_gap) && min_gap >= 0 ? Float64(min_gap) :
          throw(ArgumentError("$what: min_gap must be :auto or a finite number >= 0, got $(repr(min_gap))"))
    fgap = face_gap === :auto ? gap :
           face_gap isa Real && isfinite(face_gap) && face_gap >= 0 ? Float64(face_gap) :
           throw(ArgumentError("$what: face_gap must be :auto or a finite number >= 0, got $(repr(face_gap))"))
    inside = boundary === :inside
    for d in 1:D
        room = inside ? 2radii[1] + 2fgap : per[d] ? 2radii[1] + gap : 2radii[1]
        room < L[d] || throw(ArgumentError("$what: the largest inclusion (diameter $(2radii[1])" *
                                           (inside ? " plus 2 face_gap" : per[d] ? " plus min_gap" : "") *
                                           ") does not fit in the cell along direction $d"))
    end
    # random sequential addition saturates at φ_sat (equal sizes, many inclusions); with few
    # inclusions (n <= 2) or with free faces (inclusions may stick out) it does not apply
    eff = sum(_ballvol(Val(D), r + gap / 2) for r in radii) / V
    if (radius isa Real || radius === nothing) && length(radii) >= 3 && all(per)
        eff <= _SAT[D] || throw(ArgumentError(@sprintf("%s: %d inclusions with min_gap %.3g need an effective fraction %.3f, above the random sequential addition limit %.3f (in practice about %.2f is reachable with this gap); reduce fraction, n or min_gap",
                                                        what, length(radii), gap, eff, _SAT[D], 0.9 * _SAT[D] / (1 + gap / (2radii[end]))^D)))
    end
    obst = _obstacles(Val(D), avoid, origin, Tuple(Float64.(L)), per, radii[1] + gap)
    propose = nothing
    if cluster !== nothing
        propose, nc, Rc, share = _cluster(Val(D), cluster, rng, origin, Tuple(Float64.(L)), per, what)
        loc = share * sum(_ballvol(Val(D), r + gap / 2) for r in radii) / (nc * _ballvol(Val(D), Rc + radii[1]))
        length(radii) >= 3 && loc > _SAT[D] &&
            throw(ArgumentError(@sprintf("%s: the clusters would need a local fraction %.3f, above the random sequential addition limit %.3f; increase the cluster radius or n, or lower share",
                                         what, loc, _SAT[D])))
    end
    X = _rsa(rng, Tuple(Float64.(L)), per, radii, gap, fgap, inside, max_tries; obst, propose)
    if length(X) < length(radii)
        got = sum(_ballvol(Val(D), r) for r in radii[1:length(X)]; init=0.0) / V
        want = sum(_ballvol(Val(D), r) for r in radii) / V
        msg = @sprintf("%s: placed %d of %d inclusions (fraction %.3f of the requested %.3f) after %d attempts for the next one. Random sequential addition reaches about %.2f for equal %s with this min_gap; reduce fraction, min_gap or face_gap, or use a size distribution",
                       what, length(X), length(radii), got, want, max_tries, 0.9 * _SAT[D] / (1 + gap / (2radii[end]))^D, what == "random_spheres" ? "spheres" : what == "random_fibres" ? "fibres" : what == "particle_chains" ? "chains" : "disks")
        strict ? error(msg) : @warn(msg)
    end
    return X, radii[1:length(X)]
end

"""
    random_spheres(cell; radius, fraction, n, min_gap = :auto, face_gap = :auto,
                   boundary = :periodic, rng = Random.default_rng(), max_tries = 50_000,
                   strict = true) -> Vector{Sphere}

Non-overlapping spheres placed by random sequential addition. Give exactly two of `radius`,
`fraction` and `n`; or `radius` as a vector of radii; or `radius` as a function `rng -> r`
(size distribution) with `fraction` or `n`.

- `min_gap`: minimum distance between sphere surfaces (default `0.1 *` the smallest radius).
- `face_gap` (default `min_gap`): spheres nearly tangent to a face, edge or corner of the cell,
  or with the centre nearly on a face plane, are rejected (these give slivers and Gmsh
  periodicity failures). `face_gap = 0` switches it off.
- `boundary = :periodic`: spheres cross the faces and `generate` completes them; `:inside`:
  every sphere stays `face_gap` away from the faces (no cut spheres; leaves a depleted layer
  near the walls). Along non-periodic directions there is no wrap-around: with `:periodic` the
  spheres may cross those faces and are cut there, so the fraction inside the cell is somewhat
  below `fraction` (use `:inside` to keep them whole).
- With a radius sampler and `fraction`, radii are drawn until the next one would exceed the
  target: the packing stays at or below `fraction`. `volume_fraction(items, cell)` gives the
  value of the returned inclusions.
- `avoid`: geometries already placed (another population, `particle_chains`, items made by
  hand; a geometry, a vector, or a tuple of vectors) that the new inclusions keep `min_gap`
  away from, periodic images included. Exact for spheres, disks, capsules and fibres parallel
  to a cell axis; cylinders, ellipsoids and ellipses are replaced by the capsule that contains
  them. Mixtures: one call and one phase per population, largest first, each `fraction` of
  the whole cell: `soft = random_spheres(cell; radius = 0.05, fraction = 0.05, avoid = hard)`.
- `cluster = (n = 4, radius = 0.2)`: clustered inclusions. `n` cluster centres are drawn
  uniformly (or given as `centers = [...]`) and each inclusion is placed uniformly in the ball
  of radius `radius` around one of them; with `share = s < 1` only a fraction `s` of the
  attempts go to the clusters and the rest to the whole cell.
- `rng`: an `AbstractRNG` or an integer seed; same seed and Julia version, same packing.

Achievable fractions with equal spheres: about 0.35 (`min_gap = 0`), 0.30 (`0.1 r`), 0.26
(`0.2 r`). Above that an error reports what was placed (`strict = false`: a warning and the
partial packing). Put the result in a phase without `fuse`:
`Phase(:particles, random_spheres(cell; radius = 0.1, fraction = 0.2, rng = 1))`.
"""
function random_spheres(cell::Cell{3}; radius=nothing, fraction=nothing, n=nothing, min_gap=:auto, face_gap=:auto,
                        boundary=:periodic, avoid=nothing, cluster=nothing, rng=Random.default_rng(), max_tries=50_000, strict=true)
    per = cell.periodic
    X, r = _packing(Val(3), "random_spheres", Tuple(cell.size), per, volume(cell), radius, fraction, n,
                    min_gap, face_gap, boundary, rng, max_tries, strict; avoid, cluster, origin=Tuple(cell.origin))
    return [Sphere(cell.origin .+ x, ri) for (x, ri) in zip(X, r)]
end

"""
    random_disks(cell; radius, fraction, n, ...) -> Vector{Disk}

2D version of [`random_spheres`](@ref), same keywords. Achievable fractions with equal disks:
about 0.50 (`min_gap = 0`), 0.48 (`0.1 r`), 0.45 (`0.2 r`).
"""
function random_disks(cell::Cell{2}; radius=nothing, fraction=nothing, n=nothing, min_gap=:auto, face_gap=:auto,
                      boundary=:periodic, avoid=nothing, cluster=nothing, rng=Random.default_rng(), max_tries=50_000, strict=true)
    per = cell.periodic
    X, r = _packing(Val(2), "random_disks", Tuple(cell.size), per, volume(cell), radius, fraction, n,
                    min_gap, face_gap, boundary, rng, max_tries, strict; avoid, cluster, origin=Tuple(cell.origin))
    return [Disk(cell.origin .+ x, ri) for (x, ri) in zip(X, r)]
end

"""
    random_fibres(cell; axis = 3, radius, fraction, n, ...) -> Vector{Fibre}

Infinite fibres parallel to the cell axis `axis` (each given by its point on the mid-plane of
the cell), from a random disk packing of the cross-section (keywords as [`random_spheres`](@ref); `fraction` is the fibre volume fraction).
"""
function random_fibres(cell::Cell{3}; axis=3, radius=nothing, fraction=nothing, n=nothing, min_gap=:auto,
                       face_gap=:auto, boundary=:periodic, avoid=nothing, cluster=nothing, rng=Random.default_rng(),
                       max_tries=50_000, strict=true)
    axis isa Integer && !(axis isa Bool) && axis in 1:3 || throw(ArgumentError("random_fibres: axis must be 1, 2 or 3, got $axis"))
    t = Tuple(d for d in 1:3 if d != axis)
    L = (cell.size[t[1]], cell.size[t[2]])
    per = (cell.periodic[t[1]], cell.periodic[t[2]])
    # fibres to avoid become disks of the cross-section
    obs = Disk[]
    for g in _flat_items(avoid)
        g isa Fibre && abs(abs(g.direction[axis]) - 1) < 1e-12 ||
            throw(ArgumentError("random_fibres: avoid accepts only fibres parallel to axis $axis"))
        push!(obs, Disk((g.point[t[1]], g.point[t[2]]), g.radius))
    end
    X, r = _packing(Val(2), "random_fibres", L, per, L[1] * L[2], radius, fraction, n,
                    min_gap, face_gap, boundary, rng, max_tries, strict; avoid=obs, cluster,
                    origin=(cell.origin[t[1]], cell.origin[t[2]]))
    e = ntuple(d -> d == axis ? 1.0 : 0.0, 3)
    out = Fibre[]
    for (x, ri) in zip(X, r)
        p = ntuple(3) do d
            d == axis ? cell.origin[d] + cell.size[d] / 2 : cell.origin[d] + x[findfirst(==(d), t)]
        end
        push!(out, Fibre(p, e, ri))
    end
    return out
end

"""
    volume_fraction(items, cell)

Sum of the volumes (areas in 2D) of non-overlapping spheres, disks, ellipsoids, capsules, cylinders or fibres divided by the
cell volume, as returned by the random packings (fibres span the cell along their axis).
"""
function volume_fraction(items::AbstractVector{<:Geometry}, cell::Cell{D}) where {D}
    v = 0.0
    for it in items
        if it isa Sphere && D == 3
            v += 4π * it.radius^3 / 3
        elseif it isa Disk && D == 2
            v += π * it.radius^2
        elseif it isa Ellipsoid && D == 3
            v += 4π * prod(it.radii) / 3
        elseif it isa Ellipse && D == 2
            v += π * prod(it.radii)
        elseif it isa Capsule && D == 3
            v += π * it.radius^2 * _norm(it.p2 .- it.p1) + 4π * it.radius^3 / 3
        elseif it isa Cylinder && D == 3
            v += π * it.radius^2 * _norm(it.axis)
        elseif it isa Fibre && D == 3
            a = findfirst(x -> abs(abs(x) - 1) < 1e-12, it.direction)
            a === nothing && throw(ArgumentError("volume_fraction: only fibres parallel to a cell axis are supported"))
            v += π * it.radius^2 * cell.size[a]
        else
            throw(ArgumentError("volume_fraction(items, cell): items must be spheres or fibres (3D cell) or disks (2D cell), got $(nameof(typeof(it))) in a $(D)D cell"))
        end
    end
    return v / volume(cell)
end
