# Random packings of non-overlapping spheres, disks and aligned fibres (random sequential
# addition with gaps to the faces, edges and corners of the cell).

const _SAT = Dict(2 => 0.547, 3 => 0.384)          # RSA saturation of equal disks / spheres

_ballvol(::Val{3}, r) = 4π * r^3 / 3
_ballvol(::Val{2}, r) = π * r^2

_rng(rng::AbstractRNG) = rng
_rng(seed::Integer) = Random.Xoshiro(seed)
_rng(x) = throw(ArgumentError("rng must be an AbstractRNG or an integer seed, got $(typeof(x))"))

"""
Radii (decreasing) of the packing from the keyword combination; `V` is the volume (area in 2D)
of the region, `what` names the caller in errors.
"""
function _radii(::Val{D}, rng, V, radius, fraction, n, what) where {D}
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
    elseif radius isa Function
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
function _rsa(rng, L::NTuple{D,Float64}, per::NTuple{D,Bool}, radii, gap, fgap, inside::Bool, maxtries) where {D}
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
            x = ntuple(d -> lo[d] + rand(rng) * (hi[d] - lo[d]), D)
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
            push!(X, x); push!(R, r)
            push!(get!(bins, binof(x), Int[]), length(X))
            placed = true
            break
        end
        placed || break
    end
    return X
end

function _packing(::Val{D}, what, L, per, V, radius, fraction, n, min_gap, face_gap, boundary, rng, max_tries, strict) where {D}
    rng = _rng(rng)
    boundary in (:periodic, :inside) || throw(ArgumentError("$what: boundary must be :periodic or :inside, got $(repr(boundary))"))
    radii = _radii(Val(D), rng, V, radius, fraction, n, what)
    gap = min_gap === :auto ? 0.1 * radii[end] : Float64(min_gap)
    gap >= 0 || throw(ArgumentError("$what: min_gap must be >= 0"))
    fgap = face_gap === :auto ? gap : Float64(face_gap)
    fgap >= 0 || throw(ArgumentError("$what: face_gap must be >= 0"))
    2radii[1] + gap < minimum(L) ||
        throw(ArgumentError("$what: the largest inclusion (diameter $(2radii[1]) plus min_gap) does not fit in the cell"))
    eff = sum(_ballvol(Val(D), r + gap / 2) for r in radii) / V
    if radius isa Real || radius === nothing
        eff <= _SAT[D] || throw(ArgumentError(@sprintf("%s: %d inclusions with min_gap %.3g need an effective fraction %.3f, above the random sequential addition limit %.3f; reduce fraction, n or min_gap",
                                                        what, length(radii), gap, eff, _SAT[D])))
    end
    X = _rsa(rng, Tuple(Float64.(L)), per, radii, gap, fgap, boundary === :inside, max_tries)
    if length(X) < length(radii)
        got = sum(_ballvol(Val(D), r) for r in radii[1:length(X)]; init=0.0) / V
        want = sum(_ballvol(Val(D), r) for r in radii) / V
        msg = @sprintf("%s: placed %d of %d inclusions (fraction %.3f of the requested %.3f) after %d attempts for the next one. Random sequential addition saturates at about %.2f for equal %s (less with min_gap); reduce fraction, min_gap or face_gap, or use a size distribution",
                       what, length(X), length(radii), got, want, max_tries, _SAT[D], D == 3 ? "spheres" : "disks")
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
  near the walls). Along non-periodic directions there is no wrap-around.
- `rng`: an `AbstractRNG` or an integer seed; same seed and Julia version, same packing.

Achievable fractions with equal spheres: about 0.35 (`min_gap = 0`), 0.30 (`0.1 r`), 0.26
(`0.2 r`). Above that an error reports what was placed (`strict = false`: a warning and the
partial packing). Put the result in a phase without `fuse`:
`Phase(:particles, random_spheres(cell; radius = 0.1, fraction = 0.2, rng = 1))`.
"""
function random_spheres(cell::Cell{3}; radius=nothing, fraction=nothing, n=nothing, min_gap=:auto, face_gap=:auto,
                        boundary=:periodic, rng=Random.default_rng(), max_tries=50_000, strict=true)
    per = map(p -> p && boundary === :periodic, cell.periodic)
    X, r = _packing(Val(3), "random_spheres", Tuple(cell.size), per, volume(cell), radius, fraction, n,
                    min_gap, face_gap, boundary, rng, max_tries, strict)
    return [Sphere(cell.origin .+ x, ri) for (x, ri) in zip(X, r)]
end

"""
    random_disks(cell; radius, fraction, n, ...) -> Vector{Disk}

2D version of [`random_spheres`](@ref), same keywords. Achievable fractions with equal disks:
about 0.50 (`min_gap = 0`), 0.48 (`0.1 r`), 0.45 (`0.2 r`).
"""
function random_disks(cell::Cell{2}; radius=nothing, fraction=nothing, n=nothing, min_gap=:auto, face_gap=:auto,
                      boundary=:periodic, rng=Random.default_rng(), max_tries=50_000, strict=true)
    per = map(p -> p && boundary === :periodic, cell.periodic)
    X, r = _packing(Val(2), "random_disks", Tuple(cell.size), per, volume(cell), radius, fraction, n,
                    min_gap, face_gap, boundary, rng, max_tries, strict)
    return [Disk(cell.origin .+ x, ri) for (x, ri) in zip(X, r)]
end

"""
    random_fibres(cell; axis = 3, radius, fraction, n, ...) -> Vector{Fibre}

Infinite fibres parallel to the cell axis `axis`, from a random disk packing of the
cross-section (keywords as [`random_spheres`](@ref); `fraction` is the fibre volume fraction).
"""
function random_fibres(cell::Cell{3}; axis=3, radius=nothing, fraction=nothing, n=nothing, min_gap=:auto,
                       face_gap=:auto, boundary=:periodic, rng=Random.default_rng(), max_tries=50_000, strict=true)
    axis in 1:3 || throw(ArgumentError("random_fibres: axis must be 1, 2 or 3, got $axis"))
    t = Tuple(d for d in 1:3 if d != axis)
    L = (cell.size[t[1]], cell.size[t[2]])
    per = (cell.periodic[t[1]] && boundary === :periodic, cell.periodic[t[2]] && boundary === :periodic)
    X, r = _packing(Val(2), "random_fibres", L, per, L[1] * L[2], radius, fraction, n,
                    min_gap, face_gap, boundary, rng, max_tries, strict)
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

Sum of the volumes (areas in 2D) of non-overlapping spheres, disks or fibres divided by the
cell volume, as returned by the random packings (fibres span the cell along their axis).
"""
function volume_fraction(items::AbstractVector{<:Union{Sphere,Disk,Fibre}}, cell::Cell)
    v = 0.0
    for it in items
        v += it isa Sphere ? 4π * it.radius^3 / 3 :
             it isa Disk ? π * it.radius^2 :
             π * it.radius^2 * cell.size[argmax(abs.(it.direction))]
    end
    return v / volume(cell)
end
