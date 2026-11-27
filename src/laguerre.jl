# Dense polydisperse sphere packings and Laguerre (power) tessellations

# ---------------------------------------------------------------- dense packing

# vector from sphere i to (an image of) sphere j: the image k, or the minimum image
@inline function _pairvec(xi::NTuple{D,Float64}, xj::NTuple{D,Float64}, k, L, per, full) where {D}
    return ntuple(Val(D)) do d
        v = xj[d] - xi[d] + k[d] * L[d]
        (!full && per[d]) ? v - L[d] * round(v / L[d]) : v
    end
end
@inline _far(δ::NTuple{D,Float64}, s) where {D} = any(ntuple(d -> abs(δ[d]) >= s, Val(D)))

# Overlap relaxation: every overlapping pair is pushed apart along its centre line (the smaller
# sphere moves more), spheres are pushed inside along non-periodic directions. Returns the
# largest remaining overlap (penetration depth).
function _relax!(X::Vector{NTuple{D,Float64}}, r, L, per; tol, maxsweeps) where {D}
    N = length(X)
    rmax = maximum(r)
    # minimum image is exact when no sphere pair can touch two images of each other
    full = any(d -> per[d] && 4rmax >= L[d], 1:D)
    ks = vec(collect(Iterators.product(ntuple(d -> per[d] && full ? (-1:1) : (0:0), D)...)))
    disp = zeros(D, N)
    lim = tol * minimum(r)
    mo = Inf
    for sweep in 1:maxsweeps
        fill!(disp, 0.0)
        mo = 0.0
        @inbounds for i in 1:N-1, j in i+1:N, k in ks
            δ = _pairvec(X[i], X[j], k, L, per, full)
            s = r[i] + r[j]
            _far(δ, s) && continue
            dist = sqrt(sum(abs2, δ))
            ov = s - dist
            ov > 0 || continue
            mo = max(mo, ov)
            f = ov / max(dist, 1e-12 * s)
            wi = r[j]^D / (r[i]^D + r[j]^D)
            for d in 1:D
                disp[d, i] -= f * wi * δ[d]
                disp[d, j] += f * (1 - wi) * δ[d]
            end
        end
        @inbounds for i in 1:N, d in 1:D
            per[d] && continue
            lo, hi = r[i], L[d] - r[i]
            x = X[i][d]
            x < lo && (mo = max(mo, lo - x); disp[d, i] += lo - x)
            x > hi && (mo = max(mo, x - hi); disp[d, i] -= x - hi)
        end
        mo <= lim && return mo
        for i in 1:N
            xi = X[i]
            X[i] = ntuple(Val(D)) do d
                x = xi[d] + 0.6 * disp[d, i]
                per[d] ? x - L[d] * floor(x / L[d]) : x
            end
        end
    end
    return mo
end

# Growth from a dilute random state: the spheres are scaled by g < 1, relaxed, and grown by 1 %
# in volume per step until g = 1 or the relaxation stalls (jammed). Returns the positions and
# the last scale reached without overlap.
function _grow!(rng, R::Vector{Float64}, L::NTuple{D,Float64}, per; tol, maxsweeps=2000) where {D}
    N = length(R)
    V = prod(L)
    η = sum(_ballvol(Val(D), x) for x in R) / V
    g = min(1.0, (0.2 / η)^(1 / D))
    X = [ntuple(d -> per[d] ? rand(rng) * L[d] : g * R[i] + rand(rng) * (L[d] - 2g * R[i]), D) for i in 1:N]
    best = (copy(X), 0.0)
    while true
        mo = _relax!(X, g .* R, L, per; tol, maxsweeps)
        mo > tol * g * minimum(R) && break              # jammed
        best = (copy(X), g)
        g >= 1 && break
        g = min(1.0, g * 1.01^(1 / D))
    end
    return best
end

"""
    dense_spheres(cell; radius, fraction, n, min_gap = 0, rng = Random.default_rng(),
                  strict = true) -> Vector{Sphere} (3D) or Vector{Disk} (2D)

Non-overlapping spheres (disks in 2D) packed much denser than [`random_spheres`](@ref): the
spheres start small at random positions and grow while their overlaps are relaxed. Fractions
up to about 0.6 (3D) and 0.75 (2D, equal disks; more with a size distribution) are reached;
random sequential addition stalls at about 0.38 and 0.55. The spheres are returned largest
first (as in `random_spheres`). `radius`, `fraction` and `n` follow the rules of `random_spheres`: exactly two of
them, or a vector of radii, or a sampler `rng -> r` with `fraction` or `n`.

- Along periodic directions the spheres may cross the faces (`generate` completes them);
  along non-periodic directions they stay inside the cell.
- `min_gap`: minimum distance between sphere surfaces.
- If the packing jams before reaching the requested radii, `strict = true` throws an error
  with the fraction reached; `strict = false` warns and returns the spheres scaled down to it
  (keeping `min_gap`; an error if the smallest spheres would vanish).

The result is meant both as an inclusion phase and as the seeds of a [`laguerre`](@ref)
tessellation (`laguerre(cell; spheres = dense_spheres(...))`).
"""
function dense_spheres(cell::Cell{D}; radius=nothing, fraction=nothing, n=nothing, min_gap=0.0,
                       rng=Random.default_rng(), strict::Bool=true) where {D}
    r = _rng(rng)
    L = ntuple(d -> Float64(cell.size[d]), D)
    V = prod(L)
    R = _radii(Val(D), r, V, radius, fraction, n, "dense_spheres")
    min_gap isa Real && min_gap >= 0 || throw(ArgumentError("dense_spheres: min_gap must be a number >= 0"))
    gap = Float64(min_gap)
    per = cell.periodic
    for d in 1:D
        2maximum(R) + gap < L[d] ||
            throw(ArgumentError("dense_spheres: the largest sphere (diameter $(2maximum(R))) does not fit in the cell along direction $d"))
    end
    # effective radii: half the gap and a small margin that the relaxation may leave as overlap
    tol = 1e-3
    Re = (R .+ gap / 2) .* (1 + tol)
    X, g = _grow!(r, Re, L, per; tol)
    η = sum(_ballvol(Val(D), x) for x in R) / V
    Rout = g == 1 ? R : g .* (R .+ gap / 2) .- gap / 2     # jammed: the spheres scaled, min_gap kept
    if g < 1
        ok = all(>(0), Rout)
        ηo = ok ? sum(_ballvol(Val(D), x) for x in Rout) / V : 0.0
        ηe = sum(_ballvol(Val(D), x + gap / 2) for x in R) / V          # counting half the gap around each sphere
        msg = ok ? @sprintf("dense_spheres: the packing jammed at fraction %.3f (requested %.3f%s); reduce fraction or min_gap",
                            ηo, η, gap > 0 ? @sprintf(", %.3f counting min_gap", ηe) : "") :
              @sprintf("dense_spheres: the packing jammed and min_gap is larger than the smallest spheres can keep (requested fraction %.3f, %.3f counting min_gap); reduce fraction or min_gap",
                       η, ηe)
        (strict || !ok) && throw(ArgumentError(msg))
        @warn msg
    end
    G = D == 3 ? Sphere : Disk
    return [G(cell.origin .+ X[i], Rout[i]) for i in eachindex(R)]
end

# ---------------------------------------------------------------- Laguerre tessellation

_cellvol(P) = _poly_volume(P)

function _facearea(f::Vector{NTuple{3,Float64}})
    a = (0.0, 0.0, 0.0)
    for k in 2:length(f)-1
        a = a .+ _cross(f[k] .- f[1], f[k+1] .- f[1])
    end
    return _norm(a) / 2
end

# (neighbour, offset, area of the common face) of a raw cell
_cellfaces(P::_Poly3) = [(j, k, _facearea(f)) for (f, (j, k)) in zip(P.faces, P.labels)]
_cellfaces(P::_Poly2) = [(j, k, _norm(P.pts[m] .- P.pts[mod1(m + 1, length(P.pts))])) for (m, (j, k)) in enumerate(P.labels)]

# Weights whose cells have the target volumes v (damped Newton: the cell volumes are the
# gradient of a concave function of the weights; ∂V_i/∂w_j = −A_ij/(2 r_ij)).
function _fit_weights(X, w0, v, L, per; rtol=1e-9, maxit=100)
    N = length(X)
    w = copy(w0)
    cells(w) = try
        _voronoi_cells(X, L, per, 1e-11, w)
    catch e
        e isa _EmptyCell || rethrow()
        nothing
    end
    C = cells(w)
    C === nothing && (w = zeros(N); C = cells(w))            # start from the Voronoi cells
    C === nothing && throw(ArgumentError("laguerre: a Voronoi cell of the seeds is empty"))
    a = [_cellvol(P) for P in C]
    err = maximum(abs.(a .- v) ./ v)
    N == 1 && return w
    for _ in 1:maxit
        err < rtol && return w
        J = zeros(N, N)
        for i in 1:N, (j, k, A) in _cellfaces(C[i])
            (j == 0 || j == i) && continue
            c = A / (2 * _norm(X[j] .+ k .* L .- X[i]))
            J[i, i] += c
            J[i, j] -= c
        end
        δ = zeros(N)
        δ[2:end] = try
            J[2:end, 2:end] \ (v .- a)[2:end]
        catch e
            e isa LinearAlgebra.SingularException || rethrow()
            throw(ArgumentError("laguerre: the volumes could not be matched (disconnected cells)"))
        end
        t = 1.0
        while true
            wn = w .+ t .* δ
            Cn = cells(wn)
            if Cn !== nothing
                an = [_cellvol(P) for P in Cn]
                en = maximum(abs.(an .- v) ./ v)
                if en < (1 - t / 4) * err
                    w, C, a, err = wn, Cn, an, en
                    break
                end
            end
            t /= 2
            t < 1e-8 && throw(ArgumentError(@sprintf("laguerre: the volumes could not be matched (relative error %.2e)", err)))
        end
    end
    err < 1e3 * rtol || throw(ArgumentError(@sprintf("laguerre: the volumes could not be matched (relative error %.2e)", err)))
    return w
end

"""
    laguerre(cell; spheres = nothing, seeds = nothing, radii = nothing, volumes = nothing,
             periodic = cell.periodic) -> Tessellation

Periodic Laguerre (power, radical Voronoi) tessellation: the grain of seed `x_i` with radius
`r_i` is the set of points where `|x − x_i|² − r_i²` is smallest. Bigger radii give bigger
grains; with equal radii it is the Voronoi tessellation of the seeds. Give `spheres` (e.g. from
[`dense_spheres`](@ref): their centres and radii), or `seeds` (physical points) and `radii`
(default 0).

`volumes` sets the volume (area in 2D) of every grain exactly, adjusting the weights: a vector
that adds up to the volume of the cell, or `:spheres` for volumes proportional to the sphere
volumes (`r_i³`, 2D `r_i²`). A dense packing as seeds keeps the grains compact.

The result works with everything that takes a `Tessellation` (`grains`, `grain_phases`,
`grain_volumes`, [`grain_diameters`](@ref), `pores`, `wall_thickness`, `foam_struts`). A seed
may lie outside its own grain. `periodic` as in [`voronoi`](@ref).
"""
function laguerre(cell::Cell{D}; spheres=nothing, seeds=nothing, radii=nothing, volumes=nothing,
                  periodic=cell.periodic) where {D}
    per = periodic isa Bool ? ntuple(_ -> periodic, D) :
          (length(periodic) == D && all(x -> x isa Bool, periodic)) ? NTuple{D,Bool}(Tuple(periodic)) :
          throw(ArgumentError("laguerre: periodic must be a Bool or $D Bools"))
    L = ntuple(d -> Float64(cell.size[d]), D)
    if spheres !== nothing
        (seeds === nothing && radii === nothing) ||
            throw(ArgumentError("laguerre: give spheres, or seeds (and radii), not both"))
        G = D == 3 ? Sphere : Disk
        all(s -> s isa G, spheres) || throw(ArgumentError("laguerre: spheres must be a vector of $(G)s"))
        P = [s.center for s in spheres]
        rad = [s.radius for s in spheres]
    else
        seeds === nothing && throw(ArgumentError("laguerre: give spheres or seeds"))
        P = collect(seeds)
        rad = radii === nothing ? zeros(length(P)) :
              radii isa Real ? fill(Float64(radii), length(P)) : collect(Float64, radii)
        length(rad) == length(P) || throw(ArgumentError("laguerre: radii must give one radius per seed ($(length(P)))"))
        all(x -> isfinite(x) && x >= 0, rad) || throw(ArgumentError("laguerre: radii must be finite and >= 0"))
    end
    isempty(P) && throw(ArgumentError("laguerre: no seeds"))
    X = NTuple{D,Float64}[]
    for s in P
        p = _vec(Val(D), s, "laguerre seed") .- cell.origin
        push!(X, ntuple(D) do d
            per[d] ? p[d] - L[d] * floor(p[d] / L[d]) :
            (0 <= p[d] <= L[d] ? p[d] : throw(ArgumentError("laguerre: seed $s is outside the cell along the non-periodic direction $d")))
        end)
    end
    for i in eachindex(X), j in i+1:length(X)
        δ = ntuple(d -> (x = X[i][d] - X[j][d]; per[d] ? x - L[d] * round(x / L[d]) : x), D)
        _norm(δ) > 1e-6 * maximum(L) || throw(ArgumentError("laguerre: seeds $i and $j coincide (closer than 1e-6 of the cell)"))
    end
    w = rad .^ 2
    if volumes !== nothing
        V = prod(L)
        v = volumes === :spheres ? (any(iszero, rad) ? throw(ArgumentError("laguerre: volumes = :spheres needs positive radii")) :
                                    rad .^ D ./ sum(rad .^ D) .* V) :
            volumes isa AbstractVector ? collect(Float64, volumes) :
            throw(ArgumentError("laguerre: volumes must be a vector or :spheres"))
        length(v) == length(X) || throw(ArgumentError("laguerre: volumes must give one volume per seed ($(length(X)))"))
        all(x -> x > 0, v) || throw(ArgumentError("laguerre: volumes must be positive"))
        abs(sum(v) / V - 1) < 1e-6 ||
            throw(ArgumentError(@sprintf("laguerre: the volumes add up to %.6g, not to the volume of the cell %.6g", sum(v), V)))
        v .*= V / sum(v)                                   # exact sum (the fit matches every cell to 1e-9)
        w = _fit_weights(X, w, v, L, per)
    end
    return _tessellate(cell, X, L, per, w; what="laguerre")
end

"""
    grain_diameters(t::Tessellation) -> Vector{Float64}

Equivalent diameter of each grain: of the sphere with its volume (3D) or of the disk with its
area (2D).
"""
grain_diameters(t::Tessellation{3}) = (6 .* grain_volumes(t) ./ π) .^ (1 / 3)
grain_diameters(t::Tessellation{2}) = sqrt.(4 .* grain_volumes(t) ./ π)
