# Closed-cell foams with a pore size distribution on a Laguerre tessellation

"""
    TruncatedLogNormal

Log-normal distribution of sizes truncated to `[min, max]`, made by
[`truncated_lognormal`](@ref). Calling it with a random number generator draws a size:
`d = dist(rng)`.
"""
struct TruncatedLogNormal
    mean::Float64
    cv::Float64
    min::Float64
    max::Float64
    μ::Float64
    σ::Float64
end

Base.show(io::IO, d::TruncatedLogNormal) =
    print(io, "truncated_lognormal(", d.mean, "; cv = ", d.cv, ", min = ", d.min, ", max = ", d.max, ")")

function (d::TruncatedLogNormal)(rng)
    d.σ == 0 && return d.mean
    rng = _rng(rng)
    while true
        x = exp(d.μ + d.σ * randn(rng))
        d.min <= x <= d.max && return x
    end
end

# mean of exp(μ + σ z), z standard normal, restricted to [lo, hi]; and the probability of the window
function _tln_mean(μ, σ, lo, hi)
    # the normal density is negligible beyond ±12: integrate only where it lives
    a, b = max((log(lo) - μ) / σ, -12.0), min((log(hi) - μ) / σ, 12.0)
    a < b || return (exp(μ), 0.0)
    m = 2000
    h = (b - a) / m
    s0 = s1 = 0.0
    for k in 0:m
        z = a + k * h
        w = (k == 0 || k == m) ? 1 : (isodd(k) ? 4 : 2)
        φ = exp(-z^2 / 2)
        s0 += w * φ
        s1 += w * φ * exp(μ + σ * z)
    end
    return s1 / s0, s0 * h / 3 / sqrt(2π)
end

"""
    truncated_lognormal(mean; cv = 0.35, min, max) -> TruncatedLogNormal

Log-normal distribution of sizes (e.g. pore diameters) truncated to `[min, max]`, with the
given `mean` of the truncated distribution and `cv` the coefficient of variation of the
underlying log-normal law. `cv = 0` gives the constant `mean`. Use it as a sampler,
`d = dist(rng)` (an `AbstractRNG` or an integer seed), e.g.
`laguerre_foam(cell; diameter = truncated_lognormal(4.0; min = 2, max = 8), ...)`. In a foam the
bounds apply to the requested diameters; the measured pores are slightly smaller (see
[`pore_statistics`](@ref)).
"""
function truncated_lognormal(mean::Real; cv::Real=0.35, min::Real, max::Real)
    lo, hi, m = Float64(min), Float64(max), Float64(mean)
    isfinite(lo) && isfinite(hi) || throw(ArgumentError("truncated_lognormal: min and max must be finite"))
    0 < lo < hi || throw(ArgumentError("truncated_lognormal: need 0 < min < max"))
    lo < m < hi || throw(ArgumentError("truncated_lognormal: the mean must be between min and max"))
    cv >= 0 || throw(ArgumentError("truncated_lognormal: cv must be >= 0"))
    cv == 0 && return TruncatedLogNormal(m, 0.0, lo, hi, log(m), 0.0)
    σ = sqrt(log(1 + cv^2))
    # the truncated mean increases with μ
    a, b = log(lo) - 6σ, log(hi) + 6σ
    for _ in 1:200
        μ = (a + b) / 2
        _tln_mean(μ, σ, lo, hi)[1] < m ? (a = μ) : (b = μ)
        b - a < 1e-12 && break
    end
    μ = (a + b) / 2
    p = _tln_mean(μ, σ, lo, hi)[2]
    p > 1e-4 || throw(ArgumentError(@sprintf("truncated_lognormal: the law keeps a fraction %.2g of its mass in [min, max] (a mean very close to min or max, or a narrow window for this cv); change cv, mean or the window", p)))
    return TruncatedLogNormal(m, Float64(cv), lo, hi, μ, σ)
end

"""
    Foam{D}

Closed-cell foam made by [`laguerre_foam`](@ref): `cell`, `tessellation`, `porosity`, `wall`,
`rounding`, the target pore `diameters` (one per grain) and `scale`, the factor applied to the
requested diameters so that the pores fill the porosity of the cell. [`pores`](@ref)`(f)` gives the
pores and [`pore_statistics`](@ref)`(f)` the measured sizes.
"""
struct Foam{D}
    cell::Cell{D}
    tessellation::Tessellation{D}
    porosity::Float64
    wall::Float64
    rounding::Float64
    diameters::Vector{Float64}
    scale::Float64                  # requested diameters × scale = diameters (the pores fill the porosity)
end

Base.show(io::IO, f::Foam{D}) where {D} =
    print(io, "Foam{$D}($(length(f.diameters)) pores, porosity ", round(f.porosity; sigdigits=4), ", wall ",
          round(f.wall; sigdigits=4), ", rounding ", round(f.rounding; sigdigits=4), ")")

_check_diameters(d) = (all(x -> x > 0 && isfinite(x), d) ||
                       throw(ArgumentError("laguerre_foam: diameters must be positive and finite")); d)

_porevol(::Val{3}, d) = π * d^3 / 6
_porevol(::Val{2}, d) = π * d^2 / 4

function _draw_diameters(::Val{D}, rng, diameter, target) where {D}
    diameter isa AbstractVector && return _check_diameters(collect(Float64, diameter))
    draw = diameter isa Real ? (_ -> Float64(diameter)) : diameter
    d = Float64[]
    vol = 0.0
    while true
        x = Float64(draw(rng))
        x > 0 && isfinite(x) || throw(ArgumentError("laguerre_foam: a diameter must be positive and finite, got $x"))
        # stop at the count whose total pore volume is closest to the target
        if vol + _porevol(Val(D), x) > target
            target - vol < vol + _porevol(Val(D), x) - target || (push!(d, x))
            break
        end
        push!(d, x)
        vol += _porevol(Val(D), x)
        length(d) > 10^5 && throw(ArgumentError("laguerre_foam: more than 10^5 pores; check the diameters and the cell"))
    end
    return d
end

"""
    laguerre_foam(cell; diameter, porosity, wall = nothing, rounding = nothing,
                  packing = 0.55, rng = Random.default_rng()) -> Foam
    laguerre_foam(; n, diameter, porosity, ...)                       # cubic cell for n pores

Closed-cell foam (e.g. a metal foam) whose pores follow a size distribution. `diameter` is the
pore diameter (equivalent volume; 2D: equivalent area): a number, a vector (one per pore) or a
sampler `rng -> d` such as [`truncated_lognormal`](@ref). Diameters are drawn until the pores
fill `porosity` of the cell (with `n`, `n` diameters are drawn and the cubic cell is sized to
hold them). The pore centres are a dense packing of spheres with those relative sizes
([`dense_spheres`](@ref), fraction `packing`), and the cells of the [`laguerre`](@ref)
tessellation get volumes proportional to the pore volumes.

The walls: by default polyhedral pores (`rounding = 0`) and the wall thickness that gives the
porosity; `wall = w` keeps that wall and rounds the pores to reach the porosity
([`pore_rounding`](@ref)); `rounding = r` keeps that rounding and solves the wall. The pore
sizes drift slightly from the targets (a uniform wall removes relatively more from small
pores): [`pore_statistics`](@ref) gives both. Mesh with
`RVE(f.cell, Phase(:pores, pores(f); void = true); matrix = :walls)`.
"""
function laguerre_foam(cell::Cell{D}; diameter, porosity::Real, wall=nothing, rounding=nothing,
                       packing::Real=0.55, rng=Random.default_rng()) where {D}
    0 < porosity < 1 || throw(ArgumentError("laguerre_foam: porosity must be in (0, 1)"))
    0 < packing < (D == 3 ? 0.62 : 0.8) || throw(ArgumentError("laguerre_foam: packing must be in (0, $(D == 3 ? 0.62 : 0.8))"))
    wall === nothing || rounding === nothing || throw(ArgumentError("laguerre_foam: give wall or rounding, not both"))
    r = _rng(rng)
    V = prod(cell.size)
    d = _draw_diameters(Val(D), r, diameter, porosity * V)
    return _laguerre_foam(cell, d, Float64(porosity), wall, rounding, Float64(packing), r)
end

function laguerre_foam(; n::Integer, diameter, porosity::Real, wall=nothing, rounding=nothing,
                       packing::Real=0.55, rng=Random.default_rng())
    n >= 2 || throw(ArgumentError("laguerre_foam: n must be at least 2"))
    0 < porosity < 1 || throw(ArgumentError("laguerre_foam: porosity must be in (0, 1)"))
    0 < packing < 0.62 || throw(ArgumentError("laguerre_foam: packing must be in (0, 0.62)"))
    wall === nothing || rounding === nothing || throw(ArgumentError("laguerre_foam: give wall or rounding, not both"))
    r = _rng(rng)
    draw = diameter isa Real ? (_ -> Float64(diameter)) : diameter
    d = _check_diameters(diameter isa AbstractVector ? collect(Float64, diameter) : [Float64(draw(r)) for _ in 1:n])
    length(d) == n || throw(ArgumentError("laguerre_foam: give n diameters"))
    L = (sum(_porevol(Val(3), x) for x in d) / porosity)^(1 / 3)
    return _laguerre_foam(Cell((L, L, L)), d, Float64(porosity), wall, rounding, Float64(packing), r)
end

function _laguerre_foam(cell::Cell{D}, d, φ, wall, rounding, packing, rng) where {D}
    length(d) >= 2 || throw(ArgumentError("laguerre_foam: the cell holds fewer than 2 pores; reduce the diameters"))
    V = prod(cell.size)
    sort!(d; rev=true)                          # dense_spheres returns the spheres largest first
    pv = [_porevol(Val(D), x) for x in d]
    # sphere radii in proportion to the diameters, packed to `packing`
    # the packing only places the seeds: if it jams (few large pores), pack looser
    s = nothing
    η = packing
    for _ in 1:8
        k = (η * V / sum(_ballvol(Val(D), x / 2) for x in d))^(1 / D)
        s = try
            dense_spheres(cell; radius=k .* d ./ 2, rng)
        catch e
            e isa ArgumentError || rethrow()
            occursin("does not fit", e.msg) && throw(ArgumentError("laguerre_foam: the largest pore does not fit in the cell; enlarge the cell or reduce the diameters"))
            occursin("jammed", e.msg) || rethrow()
            nothing
        end
        s === nothing || break
        η *= 0.85
    end
    s === nothing && throw(ArgumentError("laguerre_foam: could not place the pores; check the diameters and the cell"))
    t = laguerre(cell; spheres=s, volumes=pv ./ sum(pv) .* V)
    scale = (φ * V / sum(pv))^(1 / D)
    target = d .* scale                         # the pores fill the porosity exactly
    abs(scale - 1) > 0.05 && @warn @sprintf("laguerre_foam: the requested diameters are scaled by %.3f so that the pores fill the porosity of the cell (few pores, or diameters that do not match the cell)", scale)
    if wall !== nothing
        wall isa Real && wall > 0 || throw(ArgumentError("laguerre_foam: wall must be a positive number"))
        w = Float64(wall)
        ro = try
            pore_rounding(t, 1 - φ; wall=w)
        catch e
            e isa ArgumentError || rethrow()
            throw(ArgumentError(@sprintf("laguerre_foam: porosity %.4g with wall = %.4g is out of reach (%s)", φ, w,
                                         occursin("below", e.msg) ? "the wall is too thick even for polyhedral pores" :
                                                                    "too thin even for fully rounded pores")))
        end
    else
        ro = rounding === nothing ? 0.0 : _check_rounding(rounding, "laguerre_foam")
        w = wall_thickness(t, 1 - φ; rounding=ro)
    end
    lost = count(i -> (P = _shrunk(t, i, w / 2); P === nothing || _poly_volume(P) <= 1e-9 * V), eachindex(t.loops))
    lost > 0 && @warn "laguerre_foam: $lost of $(length(d)) pores vanish with wall = $(round(w; sigdigits = 3)) (the smallest grains are thinner than the wall); use a narrower size distribution or a higher porosity"
    return Foam{D}(cell, t, φ, w, ro, target, scale)
end

pores(f::Foam) = pores(f.tessellation; wall=f.wall, rounding=f.rounding)

"""
    pore_statistics(f::Foam) -> NamedTuple

Measured pore sizes of the foam: `n` pores (`vanished`: grains thinner than the wall, left
out), `diameters` (equivalent diameter of each pore), their `mean`, `volume_weighted_mean`
(the mean a tomography sees), `std`, `min`, `max`, the `porosity` of the geometry, `wall`,
`rounding`, `faces_per_pore` (mean number of walls of a pore with other pores, periodic images
counted), `scale` and `target` (mean, min and max of the requested diameters times `scale`).
A uniform wall removes relatively more from small pores, so the measured diameters are a few
percent below the targets for the smallest pores; with few pores the drawn diameters (and
`scale`) depart from the law.
"""
function pore_statistics(f::Foam{D}) where {D}
    t = f.tessellation
    V = prod(f.cell.size)
    pv = [_pore_volume(t, i, f.wall, f.rounding) for i in eachindex(t.loops)]
    keep = pv .> 1e-9 * V
    d = D == 3 ? (6 .* pv[keep] ./ π) .^ (1 / 3) : sqrt.(4 .* pv[keep] ./ π)
    m = sum(d) / length(d)
    nf = [count(((j, _),) -> j != 0 && j != i, t.neighbours[i]) for i in eachindex(t.loops)]
    return (n=length(d), vanished=count(!, keep), diameters=d, mean=m,
            volume_weighted_mean=sum(d .^ (D + 1)) / sum(d .^ D),
            std=sqrt(sum(abs2, d .- m) / max(1, length(d) - 1)),
            min=minimum(d), max=maximum(d), porosity=sum(pv) / V, wall=f.wall, rounding=f.rounding,
            faces_per_pore=sum(nf) / length(nf), scale=f.scale,
            target=(mean=sum(f.diameters) / length(f.diameters), min=minimum(f.diameters), max=maximum(f.diameters)))
end
