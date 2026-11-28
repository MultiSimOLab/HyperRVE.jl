# Straight holes drilled through a closed-cell foam: perforated metal foams for acoustics

"""
    perforate(f::Foam; diameter, spacing = nothing, axis = 3, position = :center,
              windows = nothing, isolated = :closed, rng = Random.default_rng())
        -> (; fluid, closed, holes, windows, report)

Drill straight through-holes of `diameter` (cell units) parallel to the cell axis `axis`, on a
rectangular lattice in the plane normal to it (closed-cell foams such as Alporas absorb sound
once perforated; the hole diameter and spacing then set the flow resistivity). The holes are
[`Fibre`](@ref)s, so they are periodic along `axis` (or run between free faces) and cut the
solid; every pore a hole crosses joins the fluid.

- `spacing`: `nothing` (one hole per cell, the periodic unit cell of a perforated foam), a
  number or two numbers (the two in-plane directions in increasing order); each must divide
  the cell size (an error gives the nearest values that do). `diameter` is at most 0.8 spacing.
- `position`: `:center` (holes at the centres of the lattice cells), `:random` (the whole
  lattice shifted at random) or the in-plane coordinates `(u, v)` of one hole relative to the
  cell origin. A hole axis between 0.5 and 1.25 radii from a cell face parallel to it does not
  mesh (Gmsh), so such positions are refused (and avoided by `:random`).
- `windows`: `nothing` or the keywords of [`open_pores`](@ref) as a named tuple, e.g.
  `(open_fraction = 0.2, window = 0.4)`, to open windows in the walls as well; the
  connectivity is computed once with pores, windows and holes.
- `isolated = :closed`: the fluid is every network of pores (and windows) connected to a hole;
  the other pores are returned in `closed`. `:keep` puts every pore in `fluid`.

A pore counts as crossed when a hole goes clearly into it (to 0.95 of its radius); a pore it
only grazes stays as it is (`report.grazed_pores`). `fluid` (pores, windows and holes, which
overlap) and `closed` are meshed as with `open_pores`, e.g. `RVE(f.cell, foam_phases(a);
matrix = :aluminium)`, or `RVE(f.cell, Phase(:air, a.fluid); matrix = nothing)` for the fluid
alone; do not fuse it.

`report`: the holes, their diameter, spacing, axis, positions (absolute) and the perforation
ratio (hole area over the cell section), the pores crossed and grazed, the separate fluid
networks (a warning if more than one), whether the fluid percolates along each periodic
direction, the open and closed pores, the open porosity (pores only) and the porosity, and
the windows with their diameters and lengths.
"""
function perforate(f::Foam{3}; diameter::Real, spacing=nothing, axis::Integer=3, position=:center,
                   windows=nothing, isolated::Symbol=:closed, rng=Random.default_rng())
    axis in 1:3 || throw(ArgumentError("perforate: axis must be 1, 2 or 3"))
    diameter > 0 && isfinite(diameter) || throw(ArgumentError("perforate: diameter must be positive and finite"))
    isolated in (:closed, :keep) || throw(ArgumentError("perforate: isolated must be :closed or :keep"))
    r = _rng(rng)
    t = f.tessellation
    L = t.cell.size
    per = t.cell.periodic
    N = length(t.loops)
    ab = Tuple(d for d in 1:3 if d != axis)
    sp = _hole_spacing(spacing, L, ab)
    rh = diameter / 2
    diameter <= 0.8 * minimum(sp) + 1e-12 * maximum(L) ||
        throw(ArgumentError(@sprintf("perforate: diameter %.4g is too large for the spacing %s (at most 0.8 × spacing)", diameter, sp)))
    dm = sum(f.diameters) / length(f.diameters)
    diameter > dm && @warn @sprintf("perforate: the holes (diameter %.4g) are wider than the mean pore (%.4g)", diameter, dm)
    u0 = _hole_origin(position, sp, rh, r)
    na = (round(Int, L[ab[1]] / sp[1]), round(Int, L[ab[2]] / sp[2]))
    pts = NTuple{3,Float64}[]
    for i in 0:na[1]-1, j in 0:na[2]-1
        p = [0.0, 0.0, 0.0]
        p[ab[1]] = u0[1] + i * sp[1]
        p[ab[2]] = u0[2] + j * sp[2]
        push!(pts, Tuple(p))
    end
    H = length(pts)
    g = _foam_pores(f, "perforate")
    adj = [Tuple{Int,NTuple{3,Int}}[] for _ in 1:N+H]
    e = ntuple(d -> d == axis ? 1 : 0, 3)
    for h in 1:H
        per[axis] && push!(adj[N+h], (N + h, e))
    end
    crossed, grazed = _hole_crossings!(adj, t, g, pts, rh, axis)
    chosen = []
    if windows !== nothing
        kw = _window_options(f, windows)
        _check_windows("perforate", kw.open_fraction, kw.window, kw.min_window, kw.select)
        _, _, chosen = _choose_windows(f, g, kw.open_fraction, kw.window, kw.min_window, kw.select, r, "perforate")
        _window_edges!(adj, chosen)
        _hole_window_edges!(adj, t, chosen, pts, rh, axis)
    end
    comp, perc, _ = _networks(adj, vcat(g.alive, trues(H)))
    hc = unique(comp[N+1:N+H])
    keep = isolated === :keep ? copy(g.alive) : [g.alive[i] && comp[i] in hc for i in 1:N]
    length(hc) > 1 &&
        @warn "perforate: the holes make $(length(hc)) separate fluid networks; open windows (windows = (open_fraction = ...,)) to connect them"
    fluidperc = ntuple(d -> per[d] && any(c -> perc[c][d], hc), 3)
    wins, wd, wl, _ = _window_cylinders(t, chosen, keep)
    fluid, closed = _split_pores(g.items, keep)
    append!(fluid, wins)
    holes = [Fibre(t.cell.origin .+ p, e, rh) for p in pts]
    append!(fluid, holes)
    V = prod(L)
    report = (holes=H, diameter=Float64(diameter), spacing=sp, axis=Int(axis), positions=[t.cell.origin .+ p for p in pts],
              perforation_ratio=H * π * rh^2 / (L[ab[1]] * L[ab[2]]), crossed_pores=count(crossed),
              grazed_pores=count(grazed .& .!crossed), networks=length(hc), percolates=fluidperc,
              open_pores=count(keep), closed_pores=count(g.alive) - count(keep),
              open_porosity=sum(g.pv[keep]; init=0.0) / V, porosity=sum(g.pv) / V,
              windows=length(wins), window_diameters=wd, window_lengths=wl)
    return (fluid=fluid, closed=closed, holes=holes, windows=wins, report=report)
end

perforate(::Foam{2}; kwargs...) = throw(ArgumentError("perforate: only 3D foams"))

function _hole_spacing(spacing, L, ab)
    s = spacing === nothing ? (L[ab[1]], L[ab[2]]) :
        spacing isa Real ? (Float64(spacing), Float64(spacing)) :
        length(spacing) == 2 && all(x -> x isa Real, spacing) ? (Float64(spacing[1]), Float64(spacing[2])) :
        throw(ArgumentError("perforate: spacing must be nothing, a number or two numbers"))
    for (k, d) in enumerate(ab)
        s[k] > 0 && isfinite(s[k]) || throw(ArgumentError("perforate: spacing must be positive and finite"))
        m = L[d] / s[k]
        if round(m) < 1 || abs(m - round(m)) > 1e-6 * m
            near = unique([L[d] / max(1, floor(m)), L[d] / max(1, ceil(m))])
            throw(ArgumentError(@sprintf("perforate: spacing %.6g does not divide the cell size %.6g along direction %d; nearest: %s",
                                         s[k], L[d], d, join((@sprintf("%.6g", x) for x in near), " or "))))
        end
    end
    return s
end

# a hole axis between 0.5 and 1.25 radii from a cell face parallel to it does not mesh
_near_face(u, s, rh) = any(x -> 0.5rh < x < 1.25rh, (u, s - u))

function _hole_origin(position, sp, rh, r)
    if position === :center
        return (sp[1] / 2, sp[2] / 2)
    elseif position === :random
        for _ in 1:10_000
            u = (rand(r) * sp[1], rand(r) * sp[2])
            _near_face(u[1], sp[1], rh) || _near_face(u[2], sp[2], rh) || return u
        end
        throw(ArgumentError("perforate: no random position keeps the holes away from the cell faces"))
    elseif (position isa Tuple || position isa AbstractVector) && length(position) == 2 && all(x -> x isa Real && isfinite(x), position)
        u = (mod(Float64(position[1]), sp[1]), mod(Float64(position[2]), sp[2]))
        (_near_face(u[1], sp[1], rh) || _near_face(u[2], sp[2], rh)) &&
            throw(ArgumentError(@sprintf("perforate: a hole axis at %.4g–%.4g (0.5–1.25 radii) from a cell face does not mesh; move the holes onto the face or further from it",
                                         0.5rh, 1.25rh)))
        return u
    end
    throw(ArgumentError("perforate: position must be :center, :random or two coordinates (u, v)"))
end

function _window_options(f, windows)
    windows isa NamedTuple || throw(ArgumentError("perforate: windows must be nothing or a named tuple of open_pores keywords"))
    bad = setdiff(keys(windows), (:open_fraction, :window, :min_window, :select))
    isempty(bad) || throw(ArgumentError("perforate: unknown windows keywords $(Tuple(bad)) (open_fraction, window, min_window, select)"))
    return merge((open_fraction=0.3, window=0.4, min_window=2 * f.wall, select=:uniform), windows)
end

# lattice shifts to try: -1, 0, 1 along the periodic directions
_shifts(per) = [(a, b, c) for a in (per[1] ? (-1:1) : (0:0)), b in (per[2] ? (-1:1) : (0:0)), c in (per[3] ? (-1:1) : (0:0))]

# pores crossed (to 0.95 of the hole radius) and grazed (to 1.2 radii) by each hole, sampled on
# rings along its axis; an edge hole -> pore (with the lattice offset of the pore image) per
# crossing
function _hole_crossings!(adj, t, g, pts, rh, axis)
    L = t.cell.size
    N = length(t.loops)
    per = t.cell.periodic
    ab = Tuple(d for d in 1:3 if d != axis)
    # bounding ball of each Laguerre cell (cell-local)
    bc = Vector{NTuple{3,Float64}}(undef, N)
    bR = zeros(N)
    for i in 1:N
        g.alive[i] || continue
        V = [_vpos(t, vo) for f in t.loops[i] for vo in f]
        c = reduce(.+, V) ./ length(V)
        bc[i] = c
        bR[i] = maximum(v -> _norm(v .- c), V)
    end
    radii = (0.0, 0.25, 0.5, 0.75, 0.95, 1.05, 1.2) .* rh
    m = 16
    dz = rh / 4
    nz = max(2, ceil(Int, L[axis] / dz) + 1)
    crossed = falses(N)
    grazed = falses(N)
    for (h, p) in enumerate(pts), i in 1:N, s in _shifts(per)
        g.alive[i] || continue
        c = bc[i] .+ s .* L
        dxy = hypot(c[ab[1]] - p[ab[1]], c[ab[2]] - p[ab[2]])
        dxy <= bR[i] + 1.2rh || continue
        zlo = max(0.0, c[axis] - bR[i])
        zhi = min(L[axis], c[axis] + bR[i])
        zlo <= zhi || continue
        sL = s .* L
        best = Inf
        for z in range(0.0, L[axis]; length=nz)
            zlo - dz <= z <= zhi + dz || continue
            for ρ in radii
                ρ < best || break
                for k in 0:(ρ == 0 ? 0 : m - 1)
                    q = [p...]
                    q[axis] = z
                    q[ab[1]] += ρ * cos(2π * k / m)
                    q[ab[2]] += ρ * sin(2π * k / m)
                    if _inpore(g, Tuple(q), i, sL)
                        best = min(best, ρ)
                        break
                    end
                end
            end
        end
        isfinite(best) || continue
        grazed[i] = true
        if best <= 0.95rh
            crossed[i] = true
            push!(adj[N+h], (i, s))
            push!(adj[i], (N + h, .-s))
        end
    end
    return crossed, grazed
end

# a hole that cuts the cylinder of a window joins the two pores of the window
function _hole_window_edges!(adj, t, chosen, pts, rh, axis)
    L = t.cell.size
    N = length(t.loops)
    per = t.cell.periodic
    ab = Tuple(d for d in 1:3 if d != axis)
    for (i, j, k, c, n, rw, hi, hj) in chosen, (h, p) in enumerate(pts)
        hit = nothing
        for τ in range(-hi, hj; length=11), s in _shifts(ntuple(d -> d != axis && per[d], 3))
            q = c .+ τ .* n
            hypot(q[ab[1]] - p[ab[1]] - s[ab[1]] * L[ab[1]], q[ab[2]] - p[ab[2]] - s[ab[2]] * L[ab[2]]) < 0.95rh + rw || continue
            sa = per[axis] ? fld(q[axis], L[axis]) : 0
            hit = ntuple(d -> d == axis ? Int(sa) : s[d], 3)
            break
        end
        hit === nothing && continue
        # the hole image at `hit` touches the window of grain i (offset 0) and j (offset k)
        push!(adj[N+h], (i, .-hit))
        push!(adj[i], (N + h, hit))
        push!(adj[N+h], (j, k .- hit))
        push!(adj[j], (N + h, hit .- k))
    end
    return adj
end

"""
    foam_phases(a; air = :air, closed = :closed) -> Vector{Phase{3}}

The phases of a foam opened by [`open_pores`](@ref) or [`perforate`](@ref): the fluid as the
phase `air` and, when there are closed pores, these as the void phase `closed`, in this order
(`RVE(f.cell, foam_phases(a); matrix = :aluminium)`).
"""
function foam_phases(a; air=:air, closed=:closed)
    isempty(a.fluid) && throw(ArgumentError("foam_phases: the fluid is empty (no network percolates)"))
    ph = [Phase(air, a.fluid)]
    isempty(a.closed) || push!(ph, Phase(closed, a.closed; void=true))
    return ph
end
