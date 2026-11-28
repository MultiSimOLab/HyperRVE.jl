# Windows between the pores of a closed-cell foam: the fluid phase of an acoustic RVE

# faces between two pores, each once: grain i, neighbour j seen at the lattice offset k, the
# face polygon (cell-local coordinates of grain i) and its area
function _pore_faces(t::Tessellation{3})
    out = Tuple{Int,Int,NTuple{3,Int},Vector{NTuple{3,Float64}},Float64}[]
    z = (0, 0, 0)
    for i in eachindex(t.loops), (f, (j, k)) in zip(t.loops[i], t.neighbours[i])
        j == 0 && continue
        (i < j || (i == j && k > z)) || continue
        pts = [_vpos(t, vo) for vo in f]
        push!(out, (i, j, k, pts, _facearea(pts)))
    end
    return out
end

# largest circle inscribed in a planar convex polygon: centre (3D) and radius
function _incircle(pts::Vector{NTuple{3,Float64}})
    c0 = reduce(.+, pts) ./ length(pts)
    n = (0.0, 0.0, 0.0)
    for k in 2:length(pts)-1
        n = n .+ _cross(pts[k] .- pts[1], pts[k+1] .- pts[1])
    end
    n = n ./ _norm(n)
    u = _cross(n, abs(n[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0))
    u = u ./ _norm(u)
    v = _cross(n, u)
    q = [(_dot(p .- c0, u), _dot(p .- c0, v)) for p in pts]
    m = length(q)
    lines = Tuple{NTuple{2,Float64},Float64}[]        # inward normal a, offset b: a·x >= b
    for k in 1:m
        p, r = q[k], q[mod1(k + 1, m)]
        e = r .- p
        ne = _norm(e)
        ne > 0 || continue
        a = (-e[2], e[1]) ./ ne
        _dot(a, (0.0, 0.0) .- p) < 0 && (a = .-a)       # the centroid is inside
        push!(lines, (a, _dot(a, p)))
    end
    best, bc = 0.0, (0.0, 0.0)
    L = length(lines)
    for i in 1:L-2, j in i+1:L-1, k in j+1:L
        # a·x - ρ = b for three lines
        A = [lines[i][1][1] lines[i][1][2] -1.0; lines[j][1][1] lines[j][1][2] -1.0; lines[k][1][1] lines[k][1][2] -1.0]
        abs(LinearAlgebra.det(A)) > 1e-12 || continue
        sol = A \ [lines[i][2], lines[j][2], lines[k][2]]
        x, ρ = (sol[1], sol[2]), sol[3]
        ρ > best || continue
        all(((a, b),) -> _dot(a, x) - b >= ρ - 1e-12 * (1 + abs(ρ)), lines) || continue
        best, bc = ρ, x
    end
    return c0 .+ bc[1] .* u .+ bc[2] .* v, best, n
end

# depth along d from p at which the point enters the ball (c, R); 0 if p is inside
function _enter_depth(p, d, c, R)
    w = p .- c
    b = _dot(w, d)
    cc = _dot(w, w) - R^2
    cc <= 0 && return 0.0
    disc = b^2 - cc
    disc <= 0 && return Inf
    t = -b - sqrt(disc)
    return t >= 0 ? t : 0.0
end

"""
    open_pores(f::Foam; open_fraction = 0.3, window = 0.4, min_window = 2 f.wall,
               select = :uniform, isolated = :closed, rng = Random.default_rng())
        -> (; fluid, closed, windows, report)

Open windows in a fraction of the walls between pores (acoustic foams: a closed-cell foam
absorbs little sound; the windows and their size set the flow resistivity, the
characteristic lengths and the tortuosity). Each window is a cylinder across the wall,
normal to it, centred at the centre of the largest circle `ρ` inscribed in the wall, of radius
`window · (ρ − wall)`, and long enough for both ends to open entirely into the two pores
(it is narrowed when they do not fit; useful `window` up to about 0.6).

- `open_fraction`: the fraction of the walls between pores that get a window (over all of
  them); only walls whose window would be at least `min_window` wide are eligible (an error
  if there are too few).
- `select`: `:uniform` (each eligible wall with the same probability) or `:area`
  (probability proportional to the wall area, as membranes that break).
- `isolated = :closed`: the fluid is the largest network of pores connected through windows
  that percolates across the cell (along some periodic direction; without periodic directions,
  the largest network); every other pore stays closed and is returned in `closed` (several
  networks would make the Stokes and potential cell problems singular). `:keep` puts every
  pore in `fluid`. Windows through free faces of the cell are not made.

`fluid` (pores and windows, which overlap) is one phase: `RVE(f.cell, Phase(:air, a.fluid);
matrix = :walls)` with the solid, adding `Phase(:closed, a.closed; void = true)` when there are
closed pores (a phase needs at least one item), or `RVE(f.cell, Phase(:air, a.fluid); matrix =
nothing)` for the fluid alone (its wall is the group `interface_air_void`). Do not fuse it
(`fuse = true`): OCC fuses rounded pores and windows into surfaces whose periodic copies Gmsh
cannot match; the overlapping pieces are meshed conformally anyway. 

`report` gives the walls and those that can take a window, the windows built
(`open_fraction`: over all the walls; `open_area_fraction`: the area of the openings over the
wall area), the windows per open pore and their diameters, the connected components of all
the pores, whether the fluid percolates along each periodic direction, the open and closed
pores and the open and total porosity. If no network percolates, `fluid` is empty (with a
warning).
"""
function open_pores(f::Foam{3}; open_fraction::Real=0.3, window::Real=0.4, min_window=2 * f.wall,
                    select::Symbol=:uniform, isolated::Symbol=:closed, rng=Random.default_rng())
    0 <= open_fraction <= 1 || throw(ArgumentError("open_pores: open_fraction must be in [0, 1]"))
    0 < window < 1 || throw(ArgumentError("open_pores: window must be in (0, 1)"))
    min_window isa Real && min_window >= 0 || throw(ArgumentError("open_pores: min_window must be a number >= 0"))
    select in (:uniform, :area) || throw(ArgumentError("open_pores: select must be :uniform or :area"))
    isolated in (:closed, :keep) || throw(ArgumentError("open_pores: isolated must be :closed or :keep"))
    r = _rng(rng)
    t = f.tessellation
    L = t.cell.size
    w = f.wall
    N = length(t.loops)
    per = t.cell.periodic
    items, lost = f.rounding > 0 ? _rounded_pores(t, w, f.rounding) : _polyhedral_pores(t, w)
    lost > 0 && @warn "open_pores: $lost grains vanish with wall = $w"
    alive = [!isempty(g) for g in items]
    # the planes of every shrunk pore and, for rounded pores, its ball (cell-local)
    pplanes = Vector{Vector{Tuple{NTuple{3,Float64},Float64}}}(undef, N)
    balls = Vector{Union{Nothing,Tuple{NTuple{3,Float64},Float64}}}(nothing, N)
    for i in 1:N
        alive[i] || continue
        P = _shrunk(t, i, w / 2)
        pplanes[i] = first(_body_planes(P))
        f.rounding > 0 && (balls[i] = _pore_ball(P, f.rounding, L))
    end
    tol = 1e-9 * maximum(L)
    inpore(q, i, s) = all(((n, b),) -> _dot(n, q) <= b + _dot(n, s) + tol, pplanes[i]) &&
                      (balls[i] === nothing || _norm(q .- (balls[i][1] .+ s)) <= balls[i][2] + tol)
    # depth at which the whole end disk of a window lies inside pore i (shifted by s); Inf if none
    function depth(c, d, rw, i, s)
        for h in range(w / 2 + 0.05rw, w / 2 + 4rw; length=40)
            all(q -> inpore(q, i, s), _rim(c .+ h .* d, d, rw)) && return h
        end
        return Inf
    end
    faces = [F for F in _pore_faces(t) if alive[F[1]] && alive[F[2]]]
    isempty(faces) && throw(ArgumentError("open_pores: the foam has no walls between pores"))
    # candidate windows: centre, normal (from i to j), radius, depths into i and j; a window
    # whose end disks do not fit inside the pores is narrowed (and dropped below min_window)
    cand = []
    for (i, j, k, pts, A) in faces
        c, ρ, n = _incircle(pts)
        xj = t.seeds[j] .* L .+ k .* L
        xi = t.seeds[i] .* L
        _dot(n, xj .- xi) < 0 && (n = .-n)
        ρ > w || continue
        rw = window * (ρ - w)
        sj = k .* L
        hi = hj = Inf
        while 2rw >= min_window && rw > 0
            hi = depth(c, .-n, rw, i, (0.0, 0.0, 0.0))
            hj = depth(c, n, rw, j, sj)
            isfinite(hi) && isfinite(hj) && break
            rw *= 0.85
        end
        (2rw >= min_window && rw > 0 && isfinite(hi) && isfinite(hj)) || continue
        # no cylinder nearly tangent to a cell face (slivers, unmatched periodic faces), and
        # none through a free face of the cell
        ok = true
        for d in 1:3, x in (0.0, L[d])
            ext = ((hi + hj) / 2) * abs(n[d]) + rw * sqrt(max(0.0, 1 - n[d]^2))
            mid = c[d] + (hj - hi) / 2 * n[d]
            abs(abs(mid - x) - ext) < 0.05 * rw && (ok = false)
            !per[d] && abs(mid - x) < ext + 0.05 * rw && (ok = false)
        end
        ok || continue
        push!(cand, (i, j, k, c, n, rw, hi, hj, A))
    end
    nopen = round(Int, open_fraction * length(faces))
    nopen <= length(cand) ||
        throw(ArgumentError(@sprintf("open_pores: only %d of the %d walls between pores can take a window of at least %.3g (open fraction up to %.3f); reduce min_window or window",
                                     length(cand), length(faces), min_window, length(cand) / length(faces))))
    # :area: probability proportional to the wall area (Efraimidis–Spirakis keys log(u)/a,
    # with the areas relative to their mean so that the units do not matter)
    am = sum(c[9] for c in cand; init=0.0) / max(1, length(cand))
    keys = select === :uniform ? rand(r, length(cand)) : [log(rand(r)) / (c[9] / am) for c in cand]
    chosen = cand[sortperm(keys; rev=true)[1:nopen]]
    # connectivity with lattice offsets: components and the directions they percolate along
    adj = [Tuple{Int,NTuple{3,Int}}[] for _ in 1:N]
    for (i, j, k) in chosen
        push!(adj[i], (j, k))
        push!(adj[j], (i, .-k))
    end
    comp = zeros(Int, N)
    pos = fill((0, 0, 0), N)
    perc = NTuple{3,Bool}[]
    nc = 0
    for s in 1:N
        (alive[s] && comp[s] == 0) || continue
        nc += 1
        comp[s] = nc
        pos[s] = (0, 0, 0)
        p = [false, false, false]
        stack = [s]
        while !isempty(stack)
            a = pop!(stack)
            for (b, k) in adj[a]
                q = pos[a] .+ k
                if comp[b] == 0
                    comp[b] = nc
                    pos[b] = q
                    push!(stack, b)
                else
                    wnd = q .- pos[b]
                    for d in 1:3
                        wnd[d] != 0 && (p[d] = true)
                    end
                end
            end
        end
        push!(perc, Tuple(p))
    end
    percolating = [any(d -> per[d] && perc[c][d], 1:3) for c in 1:nc]
    # the fluid is one connected network: the largest (by pore volume) of those that percolate
    # (without periodic directions, the largest component)
    pv = [_pore_volume(t, i, w, f.rounding) for i in 1:N]
    cvol = zeros(nc)
    for i in 1:N
        comp[i] > 0 && (cvol[comp[i]] += pv[i])
    end
    cands = any(per) ? [c for c in 1:nc if percolating[c]] : collect(1:nc)
    main = isempty(cands) ? 0 : cands[argmax(cvol[cands])]
    keep = isolated === :keep ? alive : [alive[i] && comp[i] == main for i in 1:N]
    !any(keep) && @warn "open_pores: no connected network of pores percolates across the cell; raise open_fraction (all pores are returned as closed)"
    fluidperc = ntuple(d -> per[d] && main > 0 && perc[main][d], 3)
    any(keep) && any(d -> per[d] && !fluidperc[d], 1:3) &&
        @warn "open_pores: the fluid does not percolate along the periodic direction(s) $(Tuple(d for d in 1:3 if per[d] && !fluidperc[d])); raise open_fraction"
    windows = Cylinder[]
    wd = Float64[]
    nw = zeros(Int, N)
    for (i, j, k, c, n, rw, hi, hj) in chosen
        (keep[i] && keep[j]) || continue
        base = c .- hi .* n
        push!(windows, Cylinder(t.cell.origin .+ base, (hi + hj) .* n, rw))
        push!(wd, 2rw)
        nw[i] += 1
        nw[j] += 1
    end
    fluid = Geometry{3}[]
    closed = Geometry{3}[]
    for i in 1:N
        append!(keep[i] ? fluid : closed, items[i])
    end
    append!(fluid, windows)
    V = prod(L)
    areas = sum(F[5] for F in faces)
    report = (walls=length(faces), eligible=length(cand), windows=length(windows),
              open_fraction=length(windows) / length(faces), open_area_fraction=sum(π * d^2 / 4 for d in wd; init=0.0) / areas,
              windows_per_pore=sum(nw[keep]; init=0) / max(1, count(keep)), window_diameters=wd,
              components=nc, percolates=fluidperc, open_pores=count(keep), closed_pores=count(alive) - count(keep),
              open_porosity=sum(pv[keep]; init=0.0) / V, porosity=sum(pv) / V)
    return (fluid=fluid, closed=closed, windows=windows, report=report)
end

open_pores(::Foam{2}; kwargs...) = throw(ArgumentError("open_pores: only 3D foams"))

function _rim(p, n, ρ; m=16)
    u = _cross(n, abs(n[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0))
    u = u ./ _norm(u)
    v = _cross(n, u)
    return [p .+ ρ .* (cos(2π * k / m) .* u .+ sin(2π * k / m) .* v) for k in 0:m-1]
end
