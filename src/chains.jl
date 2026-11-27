# Chains of particles aligned along a cell axis (field-structured magnetorheological elastomers
# and other aligned particulate composites).

_ballvol(::Val{1}, r) = 2r

# envelope-free face bands of _rsa for one centre (used for the axial stagger and wavy chains)
function _band_ok(s, r, fgap)
    fgap == 0 && return true
    return !(abs(s - r) < fgap || (0 < s < fgap))
end

"""
    particle_chains(cell; axis = D, radius, gap, fraction = nothing, n_chains = nothing,
                    arrangement = :random, chain_gap = :auto, face_gap = :auto, fit = :gap,
                    stagger = :random, occupancy = 1.0, waviness = 0.0, modes = 2, aspect = 1.0,
                    rng = Random.default_rng(), max_tries = 50_000, strict = true)

Chains of spheres (disks in 2D) aligned along the cell axis `axis`: the microstructure of
field-structured magnetorheological elastomers. Returns the particles chain by chain, in
increasing axial coordinate; use them in a phase without `fuse`.

- `radius`: a number, a vector (one radius per chain) or a sampler `rng -> r` (one per chain);
  `gap` (> 0): minimum surface distance between consecutive particles of a chain.
- `fraction` (particle volume fraction) or `n_chains` (with `:square`/`:hex` also `(n1, n2)`);
  with a vector of radii, neither.
- `arrangement`: `:random` (random sequential addition of the chain axes in the
  cross-section), `:square` or `:hex` (in 2D `:regular`). `chain_gap` (default
  `max(gap, 0.1 r)`): minimum distance between the envelopes of two chains; `face_gap`
  (default `chain_gap`): bands to the faces, as in [`random_spheres`](@ref).
- Along a periodic axis the chains are infinite: `fit = :gap` takes as many particles as fit
  and spreads the rest of the length over the gaps (the gap only grows; a warning above 1.5×),
  `fit = :strict` requires a commensurate length. `stagger`: axial offset of each chain in
  units of the pitch (`:random`, a number in `[0, 1)`, or `:alternate` for regular
  arrangements: neighbouring chains shifted by half a pitch).
- `occupancy < 1` removes particles at random (broken chains); `waviness` (in units of the
  radius) bends the chains with `modes` smooth periodic modes (the in-chain gap never
  decreases); `aspect ≠ 1` gives spheroids (ellipses) elongated along the chain.

Reachable fractions with `gap = 0.1 r`: about 0.28 with random chains, 0.45 square, 0.52
hexagonal. `fraction` is counted before `occupancy`. Wavy chains are packed with their
amplitude as extra clearance, which leaves a thin particle-free layer next to periodic faces.
`stagger = :alternate` alternates neighbouring chains only for even counts along periodic
directions, and only along a periodic chain axis. Along free lateral faces the chains may be
cut by the face, as the inclusions of [`random_spheres`](@ref).
Mesh: the ligaments between particles of a chain have the width `gap`; `MeshOptions(gaps = 2)`
puts two elements across every thin gap and keeps the rest coarse (refining the whole phase with
`refine = gap / 2` also works but costs 10–50 times more elements).
"""
function particle_chains(cell::Cell{D}; axis=D, radius, gap, fraction=nothing, n_chains=nothing,
                         arrangement=:random, chain_gap=:auto, face_gap=:auto, fit=:gap, stagger=:random,
                         occupancy=1.0, waviness=0.0, modes=2, aspect=1.0, rng=Random.default_rng(),
                         max_tries=50_000, strict=true) where {D}
    what = "particle_chains"
    (axis isa Integer && !(axis isa Bool) && 1 <= axis <= D) || throw(ArgumentError("$what: axis must be 1…$D"))
    (gap isa Real && gap > 0) || throw(ArgumentError("$what: gap must be > 0 (for touching particles use overlapping " *
                                                     "spheres in a fused phase)"))
    fit in (:gap, :strict) || throw(ArgumentError("$what: fit must be :gap or :strict"))
    regular = D == 3 ? (:square, :hex) : (:regular, :square, :hex)
    arrangement in (:random, regular...) || throw(ArgumentError("$what: arrangement must be :random or one of $regular"))
    (occupancy isa Real && 0 < occupancy <= 1) || throw(ArgumentError("$what: occupancy must be in (0, 1]"))
    (waviness isa Real && waviness >= 0) || throw(ArgumentError("$what: waviness must be >= 0"))
    (modes isa Integer && modes >= 1) || throw(ArgumentError("$what: modes must be a positive integer"))
    (aspect isa Real && aspect > 0) || throw(ArgumentError("$what: aspect must be > 0"))
    waviness > 0 && aspect != 1 && throw(ArgumentError("$what: waviness with aspect ≠ 1 is not supported"))
    stagger === :alternate && arrangement === :random &&
        throw(ArgumentError("$what: stagger = :alternate needs a regular arrangement"))
    (stagger in (:random, :alternate) || (stagger isa Real && 0 <= stagger < 1)) ||
        throw(ArgumentError("$what: stagger must be :random, :alternate or a number in [0, 1)"))
    rng = _rng(rng)
    t = Tuple(d for d in 1:D if d != axis)
    La = Float64(cell.size[axis]); Lt = ntuple(k -> Float64(cell.size[t[k]]), D - 1)
    pera = cell.periodic[axis]; pert = ntuple(k -> cell.periodic[t[k]], D - 1)
    V = prod(cell.size)
    vpart(r) = D == 3 ? 4π * r^2 * (aspect * r) / 3 : π * r * (aspect * r)

    # --- radii per chain
    (chain_gap === :auto || (chain_gap isa Real && chain_gap >= 0)) || throw(ArgumentError("$what: chain_gap must be :auto or a number >= 0"))
    (face_gap === :auto || (face_gap isa Real && face_gap >= 0)) || throw(ArgumentError("$what: face_gap must be :auto or a number >= 0"))
    (fraction === nothing || (fraction isa Real && 0 < fraction < 1)) || throw(ArgumentError("$what: fraction must be in (0, 1)"))
    (radius isa Real || radius isa AbstractVector || !(radius isa Union{Nothing,Symbol,AbstractString})) ||
        throw(ArgumentError("$what: radius must be a number, a vector (one per chain) or a function rng -> radius"))
    if n_chains !== nothing
        ok = n_chains isa Integer ? !(n_chains isa Bool) && n_chains >= 1 :
             n_chains isa Tuple && D == 3 && arrangement !== :random && length(n_chains) == 2 &&
             all(x -> x isa Integer && !(x isa Bool) && x >= 1, n_chains)
        ok || throw(ArgumentError("$what: n_chains must be a positive integer (or (n1, n2) for a regular arrangement in 3D)"))
    end
    fit === :strict && !cell.periodic[axis] && throw(ArgumentError("$what: fit = :strict needs a periodic chain axis"))
    cg(r) = chain_gap === :auto ? max(gap, 0.1r) : Float64(chain_gap)
    fga(r) = face_gap === :auto ? cg(r) : Float64(face_gap)
    function layout(r)                         # particles per chain, pitch
        a = aspect * r
        p0 = 2a + gap
        if pera
            m = fit === :strict ? round(Int, La / p0) : floor(Int, La / p0 + 1e-12)
            m >= 1 || throw(ArgumentError("$what: one particle with its gap ($(round(p0; sigdigits=4))) does not fit along the axis"))
            fit === :strict && abs(La / p0 - round(La / p0)) > 1e-9 &&
                throw(ArgumentError("$what: fit = :strict needs the cell length along the axis to be a multiple of the pitch $p0"))
            return m, La / m
        else
            m = floor(Int, (La - 2fga(r) + gap) / p0 + 1e-12)
            m >= 1 || throw(ArgumentError("$what: no particle fits along the non-periodic axis"))
            return m, p0
        end
    end
    radii = if radius isa AbstractVector
        (fraction === nothing && n_chains === nothing) ||
            throw(ArgumentError("$what: with a vector of radii give neither fraction nor n_chains"))
        [_positive(x, "$what radius") for x in radius]
    elseif radius isa Real
        r = _positive(radius, "$what radius")
        (fraction === nothing) != (n_chains === nothing) || throw(ArgumentError("$what: give exactly one of fraction and n_chains"))
        n = if n_chains !== nothing
            n_chains isa Tuple ? prod(n_chains) : (n_chains isa Integer && n_chains >= 1 ? Int(n_chains) :
                                                  throw(ArgumentError("$what: n_chains must be a positive integer")))
        else
            (fraction isa Real && 0 < fraction < 1) || throw(ArgumentError("$what: fraction must be in (0, 1)"))
            m, _ = layout(r)
            max(1, round(Int, fraction * V / (m * vpart(r))))
        end
        fill(r, n)
    else                                        # sampler, one radius per chain
        (fraction === nothing) != (n_chains === nothing) || throw(ArgumentError("$what: give exactly one of fraction and n_chains"))
        rs = Float64[]
        if n_chains !== nothing
            for _ in 1:n_chains
                push!(rs, _positive(radius(rng), "$what sampled radius"))
            end
        else
            v = 0.0
            while true
                r = _positive(radius(rng), "$what sampled radius")
                m, _ = layout(r)
                v + m * vpart(r) > fraction * V && break
                push!(rs, r); v += m * vpart(r)
                length(rs) > 10^5 && throw(ArgumentError("$what: too many chains; check the sampler"))
            end
            isempty(rs) && throw(ArgumentError("$what: one chain already exceeds the requested fraction"))
        end
        rs
    end
    A = waviness                                  # amplitude in units of the radius

    # --- lateral positions of the chain axes (cell-local, transverse coordinates)
    X = NTuple{D - 1,Float64}[]
    R = Float64[]
    par = Int[]                                   # parity of the chain in a regular arrangement
    if arrangement === :random
        # the chain axes are packed with the particle radius; a wavy chain moves up to A r off its
        # axis, so the gaps and the face bands grow by that amount (2 A r between two chains)
        rmin, rmax = extrema(radii)
        Xs, Rs = _packing(Val(D - 1), what, Lt, pert, prod(Lt), radii, nothing, nothing, cg(rmin) + 2A * rmax,
                          fga(rmin) + A * rmax, :periodic, rng, max_tries, strict)
        X, R = Xs, Rs
    else
        n = length(radii)
        nn = if D == 2
            (n,)
        elseif n_chains isa Tuple
            (Int(n_chains[1]), Int(n_chains[2]))
        else
            n1 = max(1, round(Int, sqrt(n * Lt[1] / Lt[2])))
            n2 = max(1, round(Int, n / n1))
            n1 * n2 == n || throw(ArgumentError("$what: $n chains do not form a regular $(arrangement) arrangement of this cross-section; use n_chains = ($n1, $n2) or another number"))
            (n1, n2)
        end
        hex = arrangement === :hex && D == 3
        for j in 1:(D == 3 ? nn[2] : 1), i in 1:nn[1]
            x1 = (i - 0.5 + (hex && iseven(j) ? 0.5 : 0.0)) * Lt[1] / nn[1]
            push!(X, D == 3 ? (pert[1] ? mod(x1, Lt[1]) : x1, (j - 0.5) * Lt[2] / nn[2]) : (x1,))
            push!(par, mod(i + j, 2))
        end
        rr = radii[1]
        all(==(rr), radii) || throw(ArgumentError("$what: regular arrangements need equal radii"))
        R = fill(rr, length(X))
        if hex
            ratio = (Lt[2] / nn[2]) / (Lt[1] / nn[1])
            abs(ratio - sqrt(3) / 2) > 0.01 * sqrt(3) / 2 &&
                @warn "$what: the cross-section does not give a hexagonal arrangement (row spacing / chain spacing = $(round(ratio; digits=3)), hexagonal 0.866)"
        end
        D == 3 && hex && pert[2] && isodd(nn[2]) &&
            @warn "$what: an odd number of rows along a periodic direction cannot be hexagonal (the first and last rows are neighbours)"
        stagger === :alternate && any(k -> pert[k] && isodd(nn[k]), 1:D-1) &&
            @warn "$what: stagger = :alternate with an odd number of chains along a periodic direction: some neighbours share the offset"
        stagger === :alternate && !pera && @warn "$what: stagger = :alternate is ignored along a non-periodic chain axis"
        # every chain inside the free faces; spacing to the others and to its own images
        for xa in X, k in 1:D-1
            pert[k] && continue
            (xa[k] - rr * (1 + A) >= fga(rr) - 1e-12 && xa[k] + rr * (1 + A) <= Lt[k] - fga(rr) + 1e-12) ||
                throw(ArgumentError("$what: a chain of the regular arrangement crosses a free face; use fewer chains or a periodic cell"))
        end
        for k in 1:D-1
            pert[k] && Lt[k] < 2 * rr * (1 + A) + cg(rr) - 1e-12 &&
                throw(ArgumentError("$what: a chain is closer than chain_gap to its own periodic image"))
        end
        # spacing check (minimum image)
        for a in eachindex(X), b in a+1:length(X)
            δ = ntuple(k -> (x = X[a][k] - X[b][k]; pert[k] ? x - Lt[k] * round(x / Lt[k]) : x), D - 1)
            _norm(δ) >= 2 * rr * (1 + A) + cg(rr) - 1e-12 ||
                throw(ArgumentError("$what: the chains of the regular arrangement are closer than 2r + chain_gap; use fewer chains"))
        end
    end

    # --- particles of each chain
    items = Geometry{D}[]
    out_chains = Vector{Vector{NTuple{D,Float64}}}()
    out_r = Float64[]
    base_alt = nothing
    gap_warned = false
    band_warned = false
    for (c, (x, r)) in enumerate(zip(X, R))
        m, p = layout(r)
        a = aspect * r
        if pera && !gap_warned && p - 2a > 1.5gap
            @warn "$what: to close the chains along the periodic axis the gap grows to $(round(p - 2a; sigdigits=3)) (asked $gap)"
            gap_warned = true
        end
        fa = fga(r)
        # the axial band is narrower by default: a chain must fit its gap between two bands
        fx = face_gap === :auto ? min(fa, a / 4) : fa
        z0_ok(z0) = (s = min(mod(z0, p), p - mod(z0, p)); _band_ok(s, a, fx) || s == 0)
        # every particle against the face, edge and corner bands (periodic directions)
        function particle_ok(q)
            sd = ntuple(d -> min(q[d], cell.size[d] - q[d]), D)
            for d in 1:D
                cell.periodic[d] || continue
                rd, f = d == axis ? (a, fx) : (r, fa)
                sd[d] < 1e-12 * cell.size[d] && continue          # a clean cut through the centre
                (abs(sd[d] - rd) < f || sd[d] < f) && return false
            end
            aspect == 1 || return true
            P = [d for d in 1:D if cell.periodic[d] && sd[d] >= 1e-12 * cell.size[d]]
            f = min(fa, fx)
            for i in eachindex(P), j in i+1:length(P)
                abs(hypot(sd[P[i]], sd[P[j]]) - r) < f && return false
            end
            length(P) == 3 && abs(_norm(sd) - r) < f && return false
            return true
        end
        retry = pera ? stagger === :random || A > 0 : A > 0 || stagger === :random
        local pts
        found = false
        for attempt in 1:(retry ? max_tries : 1)
            z0 = if pera
                if stagger isa Real
                    stagger * p
                elseif stagger === :alternate
                    if base_alt === nothing
                        for _ in 1:max_tries
                            b = rand(rng) * p
                            z0_ok(b) && z0_ok(b + p / 2) && (base_alt = b; break)
                        end
                        base_alt === nothing &&
                            throw(ArgumentError("$what: no axial offset of the alternating chains satisfies the face bands; reduce face_gap or use stagger = 0.5"))
                    end
                    base_alt + (par[c] == 1 ? p / 2 : 0.0)
                else
                    z = nothing
                    for _ in 1:max_tries
                        zz = rand(rng) * p
                        z0_ok(zz) && (z = zz; break)
                    end
                    z === nothing && throw(ArgumentError("$what: no axial offset satisfies the face bands; reduce face_gap"))
                    z
                end
            else
                slack = La - 2fa - (m * 2a + (m - 1) * gap)
                u = stagger isa Real ? stagger : stagger === :random ? rand(rng) : 0.5
                fa + a + u * slack
            end
            # waviness: smooth periodic lateral offset, max |w| <= A r
            dirs = [D == 3 ? (θ = 2π * rand(rng); (cos(θ), sin(θ))) : (rand(rng) < 0.5 ? (-1.0,) : (1.0,)) for _ in 1:modes]
            ph = [2π * rand(rng) for _ in 1:modes]
            cw = [rand(rng) for _ in 1:modes]; cw ./= sum(cw)
            w(z) = A == 0 ? ntuple(_ -> 0.0, D - 1) :
                   reduce((u, v) -> u .+ v, (A * r * cw[k] * cos(2π * k * z / La + ph[k]) .* dirs[k] for k in 1:modes))
            pts = NTuple{D,Float64}[]
            for i in 0:m-1
                z = z0 + i * p
                pera && (z = mod(z, La))
                xt = x .+ w(z)
                xt = ntuple(k -> pert[k] ? mod(xt[k], Lt[k]) : xt[k], D - 1)
                push!(pts, ntuple(d -> d == axis ? z : xt[findfirst(==(d), t)], D))
            end
            if all(particle_ok, pts)
                found = true
                break
            end
        end
        if !found
            retry && throw(ArgumentError("$what: could not place chain $c inside the face, edge and corner bands; " *
                                         "reduce face_gap or waviness"))
            band_warned || @warn "$what: some particles are nearly tangent to a face, edge or corner of the cell " *
                                 "(stagger = $(repr(stagger))); expect poor elements there"
            band_warned = true
        end
        sort!(pts; by=q -> q[axis])
        push!(out_chains, pts); push!(out_r, r)
    end
    # occupancy
    keep = [[occupancy == 1 || rand(rng) < occupancy for _ in ch] for ch in out_chains]
    any(any, keep) || (keep[1][1] = true)
    for (ch, k, r) in zip(out_chains, keep, out_r), (q, kk) in zip(ch, k)
        kk || continue
        cpos = cell.origin .+ q
        if aspect == 1
            push!(items, D == 3 ? Sphere(cpos, r) : Disk(cpos, r))
        else
            radv = ntuple(d -> d == axis ? aspect * r : r, D)
            push!(items, D == 3 ? Ellipsoid(cpos, radv) : Ellipse(cpos, radv))
        end
    end
    return D == 3 ? (aspect == 1 ? Sphere[i for i in items] : Ellipsoid[i for i in items]) :
                    (aspect == 1 ? Disk[i for i in items] : Ellipse[i for i in items])
end
