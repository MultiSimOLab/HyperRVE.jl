# Triply periodic minimal surfaces (experimental): level-set solids (network f ≤ c, sheet |f| ≤ t)
# triangulated by periodic marching tetrahedra on a grid of the unit cube, with snapping of the
# grid nodes, edge flips and smoothing projected onto the surface. The facets of opposite faces of
# the cell coincide by construction. Kernel only (pure Julia): densities and levels.

const _TPMS_F = Dict{Symbol,Function}(
    :gyroid => (x, y, z) -> sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x),
    :schwarz_p => (x, y, z) -> cos(x) + cos(y) + cos(z),
    :diamond => (x, y, z) -> sin(x) * sin(y) * sin(z) + sin(x) * cos(y) * cos(z) + cos(x) * sin(y) * cos(z) + cos(x) * cos(y) * sin(z),
    :iwp => (x, y, z) -> 2 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x)) - (cos(2x) + cos(2y) + cos(2z)),
    :neovius => (x, y, z) -> 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z))

# F < 0 is the solid. kind = :network (f < level) or :sheet (|f| < level)
function _tp_levelfun(name, kind, level, n, shift=0.0)
    f = _TPMS_F[name]
    w = 2π * n
    kind == :network ? (p -> f(w * (p[1] + shift), w * (p[2] + shift), w * (p[3] + shift)) - level) :
                       (p -> abs(f(w * (p[1] + shift), w * (p[2] + shift), w * (p[3] + shift))) - level)
end

const _TP_IFACE = 0x40
const _TP_SNAP_MAXPLANES = 1   # nodes on cell edges are not snapped (pinches the surface there)
_tp_bit(d, s) = UInt8(1) << (2 * (d - 1) + s - 1)

mutable struct _TPSurf
    P::Vector{NTuple{3,Float64}}
    feat::Vector{UInt8}
    key::Vector{NTuple{2,Int}}          # grid edge (ga, gb), ga == gb for grid nodes
    tri::Vector{NTuple{3,Int}}
    lab::Vector{Int}                    # 0 interface (outward of solid); f solid cap on face f=1..6; 6+f complement cap
    N::Int
end

function _tp_march(F, N::Int; root=true, snap=0.0)
    M = N + 1
    gid(i, j, k) = 1 + i + M * (j + M * k)
    gpos(g) = (g -= 1; (g % M, (g ÷ M) % M, g ÷ (M * M)))
    V = Array{Float64}(undef, M, M, M)
    for k in 0:N, j in 0:N, i in 0:N
        v = F((mod(i, N) / N, mod(j, N) / N, mod(k, N) / N))
        V[i+1, j+1, k+1] = abs(v) < 1e-9 ? 1e-9 : v
    end
    val(g) = V[g]
    xyz(g) = gpos(g) ./ N
    featg(g) = (c = gpos(g); r = 0x00;
                for d in 1:3
                    c[d] == 0 && (r |= _tp_bit(d, 1)); c[d] == N && (r |= _tp_bit(d, 2))
                end; r)
    S = _TPSurf(NTuple{3,Float64}[], UInt8[], NTuple{2,Int}[], NTuple{3,Int}[], Int[], N)
    ids = Dict{NTuple{2,Int},Int}()
    # snapping: a grid node with a cut point closer than snap (fraction of the edge) along an edge of
    # its own feature (cell face / cell edge) is moved to that cut point; every cut point of its edges
    # collapses there. Decided in an order that is the same for periodic copies.
    snapto = zeros(Int, M^3)
    if snap > 0
        # decided per periodic class (torus node), so that the copies of a node on opposite faces agree
        nbrs = ((1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1))
        tid(c) = 1 + mod(c[1], N) + N * (mod(c[2], N) + N * mod(c[3], N))
        tval(c) = V[mod(c[1], N)+1, mod(c[2], N)+1, mod(c[3], N)+1]
        cand = Tuple{Int,Float64,Int,NTuple{3,Int},NTuple{3,Int}}[]
        for k in 0:N-1, j in 0:N-1, i in 0:N-1
            c0 = (i, j, k); a = gid(i, j, k); fa = featg(a); count_ones(fa) >= _TP_SNAP_MAXPLANES + 1 && continue
            best = (Inf, (0, 0, 0))
            for o in nbrs, sg in (1, -1)
                c = c0 .+ sg .* o
                all(x -> 0 <= x <= N, c) || continue
                b = gid(c...)
                ((val(a) < 0) != (val(b) < 0) && (featg(b) & fa) == fa) || continue
                sp = abs(val(a)) / (abs(val(a)) + abs(val(b)))
                sp < best[1] && (best = (sp, sg .* o))
            end
            best[1] < snap && push!(cand, (-count_ones(fa), best[1], tid(c0), c0, best[2]))
        end
        tsnap = Dict{Int,NTuple{3,Int}}()
        for (_, _, t, c0, o) in sort!(cand)
            blocked = false
            for q in nbrs, sg in (1, -1)
                c = c0 .+ sg .* q
                (tval(c) < 0) != (tval(c0) < 0) && haskey(tsnap, tid(c)) && (blocked = true)
            end
            blocked || (tsnap[t] = o)
        end
        for k in 0:N, j in 0:N, i in 0:N
            o = get(tsnap, tid((i, j, k)), nothing)
            o === nothing && continue
            snapto[gid(i, j, k)] = gid(((i, j, k) .+ o)...)
        end
    end
    function node(a, b)
        a > b && ((a, b) = (b, a))
        a != b && snapto[a] != 0 && return node(a, a)
        a != b && snapto[b] != 0 && return node(b, b)
        get!(ids, (a, b)) do
            if a == b && snapto[a] == 0
                push!(S.P, xyz(a)); push!(S.feat, featg(a))
            else
                fsn = a == b ? (featg(a) | _TP_IFACE) : ((featg(a) & featg(b)) | _TP_IFACE)
                ka = a
                a == b && (b = snapto[a])
                pa, pb = xyz(a), xyz(b); fa, fb = val(a), val(b)
                s = fa / (fa - fb)
                if root   # true root of F on the grid edge (bisection; the edge is periodic-consistent)
                    lo, hi = 0.0, 1.0
                    for _ in 1:40
                        m = (lo + hi) / 2
                        fm = F(pa .+ m .* (pb .- pa))
                        (fm < 0) == (fa < 0) ? (lo = m) : (hi = m)
                    end
                    s = (lo + hi) / 2
                end
                # cut points on cell edges cannot be smoothed: keep them away from the grid nodes
                s = count_ones(fsn & 0x3f) >= 2 ? clamp(s, 0.15, 0.85) : clamp(s, 1e-3, 1 - 1e-3)
                push!(S.P, pa .+ s .* (pb .- pa)); push!(S.feat, fsn)
                a = ka; fsn & 0x3f == featg(a) && snapto[a] != 0 && (b = a)
            end
            push!(S.key, (a, b))
            length(S.P)
        end
    end
    perms = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
    e(d) = ntuple(i -> i == d ? 1 : 0, 3)
    sub(a, b) = a .- b
    cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
    dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
    function emit(t, lab, dir)      # orient so that normal·dir > 0
        (t[1] == t[2] || t[2] == t[3] || t[1] == t[3]) && return
        nrm = cross(sub(S.P[t[2]], S.P[t[1]]), sub(S.P[t[3]], S.P[t[1]]))
        push!(S.tri, dot(nrm, dir) >= 0 ? t : (t[1], t[3], t[2])); push!(S.lab, lab)
    end
    for k in 0:N-1, j in 0:N-1, i in 0:N-1, p in perms
        c0 = (i, j, k); c1 = c0 .+ e(p[1]); c2 = c1 .+ e(p[2]); c3 = c2 .+ e(p[3])
        g = (gid(c0...), gid(c1...), gid(c2...), gid(c3...))
        ins = [x for x in g if val(x) < 0]; out = [x for x in g if val(x) >= 0]
        (isempty(ins) || isempty(out)) && continue
        dir = sub(xyz(out[1]), xyz(ins[1]))
        if length(ins) == 1
            emit((node(ins[1], out[1]), node(ins[1], out[2]), node(ins[1], out[3])), 0, dir)
        elseif length(ins) == 3
            emit((node(ins[1], out[1]), node(ins[2], out[1]), node(ins[3], out[1])), 0, dir)
        else
            a, b = ins; c, d = out
            q = (node(a, c), node(a, d), node(b, d), node(b, c))
            emit((q[1], q[2], q[3]), 0, dir); emit((q[1], q[3], q[4]), 0, dir)
        end
    end
    # caps: grid triangles of the six faces, clipped by the sign
    for d in 1:3, s in 1:2
        a, b = d == 1 ? (2, 3) : d == 2 ? (1, 3) : (1, 2)
        dir = s == 1 ? .-e(d) : e(d)
        face = 2 * (d - 1) + s
        for v in 0:N-1, u in 0:N-1
            c00 = ntuple(i -> i == d ? (s == 1 ? 0 : N) : i == a ? u : v, 3)
            c10 = c00 .+ e(a); c01 = c00 .+ e(b); c11 = c10 .+ e(b)
            for T in ((c00, c10, c11), (c00, c01, c11))
                g = (gid(T[1]...), gid(T[2]...), gid(T[3]...))
                for side in 0:1
                    poly = Int[]
                    for m in 1:3
                        x, y = g[m], g[mod1(m + 1, 3)]
                        inx = (val(x) < 0) == (side == 0); iny = (val(y) < 0) == (side == 0)
                        inx && push!(poly, node(x, x))
                        inx != iny && push!(poly, node(x, y))
                    end
                    for m in 2:length(poly)-1
                        emit((poly[1], poly[m], poly[m+1]), face + 6 * side, dir)
                    end
                end
            end
        end
    end
    return S, ids
end

function _tp_signed_volume(S::_TPSurf, labs)
    v = 0.0
    for (t, l) in zip(S.tri, S.lab)
        l in labs || continue
        a, b, c = S.P[t[1]], S.P[t[2]], S.P[t[3]]
        v += a[1] * (b[2] * c[3] - b[3] * c[2]) + a[2] * (b[3] * c[1] - b[1] * c[3]) + a[3] * (b[1] * c[2] - b[2] * c[1])
    end
    return v / 6
end
_tp_solid_volume(S) = _tp_signed_volume(S, 0:6)
function _tp_area(S::_TPSurf, labs=0:0)
    A = 0.0
    for (t, l) in zip(S.tri, S.lab)
        l in labs || continue
        a, b, c = S.P[t[1]], S.P[t[2]], S.P[t[3]]
        u = b .- a; w = c .- a
        A += sqrt((u[2] * w[3] - u[3] * w[2])^2 + (u[3] * w[1] - u[1] * w[3])^2 + (u[1] * w[2] - u[2] * w[1])^2) / 2
    end
    A
end

function _tp_edges(S::_TPSurf)
    ef = Dict{NTuple{2,Int},UInt8}()
    for (t, l) in zip(S.tri, S.lab), m in 1:3
        a, b = minmax(t[m], t[mod1(m + 1, 3)])
        f = (S.feat[a] & S.feat[b] & 0x3f) | (l == 0 ? _TP_IFACE : 0x00)
        ef[(a, b)] = get(ef, (a, b), 0x00) | f
    end
    ef
end

# Laplacian smoothing restricted to the feature of each node + projection onto F = 0; periodic copies stay identical.
function _tp_smooth!(S::_TPSurf, F, ids; iters=20, λ=0.6)
    N = S.N; M = N + 1
    ef = _tp_edges(S)
    nb = [Int[] for _ in S.P]
    for ((a, b), f) in ef
        (f & S.feat[a]) == S.feat[a] && push!(nb[a], b)
        (f & S.feat[b]) == S.feat[b] && push!(nb[b], a)
    end
    foreach(sort!, nb)
    free = [ntuple(d -> (S.feat[p] & (_tp_bit(d, 1) | _tp_bit(d, 2))) == 0, 3) for p in eachindex(S.P)]
    # periodic master of each node: same grid edge with indices N -> 0
    gpos(g) = (g -= 1; (g % M, (g ÷ M) % M, g ÷ (M * M)))
    gid(c) = 1 + c[1] + M * (c[2] + M * c[3])
    function canon(k)
        ca, cb = gpos(k[1]), gpos(k[2])
        sh = ntuple(d -> (ca[d] == N && cb[d] == N) ? N : 0, 3)
        a, b = gid(ca .- sh), gid(cb .- sh)
        (min(a, b), max(a, b)), sh ./ N
    end
    master = [(q = canon(S.key[p]); (ids[q[1]], q[2])) for p in eachindex(S.P)]
    δ = 1e-6
    grad(p) = ntuple(d -> (F(ntuple(i -> p[i] + (i == d ? δ : 0.0), 3)) - F(ntuple(i -> p[i] - (i == d ? δ : 0.0), 3))) / 2δ, 3)
    for it in 1:iters
        newP = copy(S.P)
        for p in eachindex(S.P)
            master[p][1] == p || continue
            isempty(nb[p]) && continue
            c = (0.0, 0.0, 0.0)
            for q in nb[p]
                c = c .+ S.P[q]
            end
            c = c ./ length(nb[p])
            x = S.P[p] .+ λ .* (c .- S.P[p]) .* free[p]
            if S.feat[p] & _TP_IFACE != 0
                for _ in 1:3
                    g = grad(x) .* free[p]
                    g2 = sum(abs2, g)
                    g2 > 1e-12 || break
                    x = x .- (F(x) / g2) .* g
                end
            end
            newP[p] = x
        end
        for p in eachindex(S.P)
            m, sh = master[p]
            m == p || (newP[p] = newP[m] .+ sh)
        end
        S.P = newP
    end
    return S
end

# each edge of the closed shell of a phase must be in exactly two triangles
function _tp_watertight(S::_TPSurf, labs)
    cnt = Dict{NTuple{2,Int},Int}()
    for (t, l) in zip(S.tri, S.lab), m in 1:3
        l in labs || continue
        e = minmax(t[m], t[mod1(m + 1, 3)])
        cnt[e] = get(cnt, e, 0) + 1
    end
    all(==(2), values(cnt))
end

# connected components (shells) of the triangles with labels in labs: comp index per triangle (0 if not in labs)
function _tp_shells(S::_TPSurf, labs)
    parent = collect(1:length(S.P))
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for (t, l) in zip(S.tri, S.lab)
        l in labs || continue
        r = find(t[1]); parent[find(t[2])] = r; parent[find(t[3])] = r
    end
    roots = Dict{Int,Int}()
    comp = zeros(Int, length(S.tri))
    for (i, (t, l)) in enumerate(zip(S.tri, S.lab))
        l in labs || continue
        comp[i] = get!(roots, find(t[1]), length(roots) + 1)
    end
    comp, length(roots)
end

_tp_minang(a, b, c) = (l = (sqrt(sum(abs2, b .- a)), sqrt(sum(abs2, c .- b)), sqrt(sum(abs2, a .- c)));
    minimum(acosd(clamp((l[m]^2 + l[mod1(m + 1, 3)]^2 - l[mod1(m + 2, 3)]^2) / (2 * l[m] * l[mod1(m + 1, 3)] + 1e-300), -1, 1)) for m in 1:3))
_tp_nrm(a, b, c) = (u = b .- a; w = c .- a; (u[2] * w[3] - u[3] * w[2], u[3] * w[1] - u[1] * w[3], u[1] * w[2] - u[2] * w[1]))

# edge flips inside the interface (never on feature curves, never in the caps: the caps of opposite
# faces must stay identical)
function _tp_flip!(S::_TPSurf; mincos=0.5)
    et = Dict{NTuple{2,Int},Vector{Int}}()
    for (i, t) in enumerate(S.tri), m in 1:3
        push!(get!(et, minmax(t[m], t[mod1(m + 1, 3)]), Int[]), i)
    end
    touched = falses(length(S.tri)); nf = 0
    for e in sort!(collect(keys(et)))
        ts = et[e]
        length(ts) == 2 || continue
        i, j = ts
        (S.lab[i] == 0 && S.lab[j] == 0) || continue
        (touched[i] || touched[j]) && continue
        a, b = e
        ti, tj = S.tri[i], S.tri[j]
        c = only(x for x in ti if x != a && x != b); d = only(x for x in tj if x != a && x != b)
        haskey(et, minmax(c, d)) && continue
        # orientation of ti: a->b or b->a
        k = findfirst(==(a), ti); fwd = ti[mod1(k + 1, 3)] == b
        n1, n2 = fwd ? ((a, d, c), (d, b, c)) : ((a, c, d), (c, b, d))
        P = S.P
        old = min(_tp_minang(P[ti[1]], P[ti[2]], P[ti[3]]), _tp_minang(P[tj[1]], P[tj[2]], P[tj[3]]))
        new = min(_tp_minang(P[n1[1]], P[n1[2]], P[n1[3]]), _tp_minang(P[n2[1]], P[n2[2]], P[n2[3]]))
        new > old + 1.0 || continue
        no = _tp_nrm(P[ti[1]], P[ti[2]], P[ti[3]]); nn1 = _tp_nrm(P[n1[1]], P[n1[2]], P[n1[3]]); nn2 = _tp_nrm(P[n2[1]], P[n2[2]], P[n2[3]])
        cs(u, v) = sum(u .* v) / sqrt(sum(abs2, u) * sum(abs2, v) + 1e-300)
        (cs(no, nn1) > mincos && cs(no, nn2) > mincos && cs(nn1, nn2) > mincos) || continue
        S.tri[i] = n1; S.tri[j] = n2; touched[i] = touched[j] = true; nf += 1
    end
    nf
end

# ---------------------------------------------------------------- surfaces, levels and densities

const _TPMS_KINDS = (:gyroid, :schwarz_p, :diamond, :iwp, :neovius)

function _tpms_args(kind, type, periods, resolution, what)
    kind in _TPMS_KINDS || throw(ArgumentError("$what: kind must be one of $(join(_TPMS_KINDS, ", ")), got $(repr(kind))"))
    type in (:sheet, :network) || throw(ArgumentError("$what: type must be :sheet or :network, got $(repr(type))"))
    periods isa Integer && !(periods isa Bool) && periods >= 1 || throw(ArgumentError("$what: periods must be a positive integer"))
    res = resolution === nothing ? 24 : resolution
    res isa Integer && !(res isa Bool) && res >= 8 || throw(ArgumentError("$what: resolution must be an integer >= 8 (grid cubes per period)"))
    return Int(periods), Int(res)
end

# range of f over a period (sampled), for the admissible levels
function _tpms_range(kind)
    f = _TPMS_F[kind]
    r = range(0, 2π; length=49)
    vals = [f(x, y, z) for x in r, y in r, z in r]
    return minimum(vals), maximum(vals)
end

# marched, smoothed and flipped surface of the solid in the unit cube (periodic in the three
# directions); shift translates the level-set function (a fraction of the cell)
function _tpms_surface(kind, type, level, n, N; snap=0.25, shift=0.0)
    F = _tp_levelfun(kind, type, level, n, shift)
    S, ids = _tp_march(F, N; snap=snap)
    for _ in 1:3
        _tp_smooth!(S, F, ids; iters=5); _tp_flip!(S)
    end
    _tp_smooth!(S, F, ids; iters=5)
    return S, F
end

# defects that the watertightness check does not see: vertices with more than one fan of triangles
# (pinched surface) and folded triangles (interface against the gradient, caps against the face)
function _tp_defects(S, F)
    bad = 0
    δ = 1e-7
    for (t, l) in zip(S.tri, S.lab)
        a, b, c = S.P[t[1]], S.P[t[2]], S.P[t[3]]
        nrm = _tp_nrm(a, b, c)
        if l == 0
            m = (a .+ b .+ c) ./ 3
            g = ntuple(d -> (F(ntuple(i -> m[i] + (i == d ? δ : 0.0), 3)) - F(ntuple(i -> m[i] - (i == d ? δ : 0.0), 3))) / 2δ, 3)
            sum(nrm .* g) > 0 || (bad += 1)
        else
            f = mod1(l, 6); d = (f + 1) ÷ 2
            nrm[d] * (isodd(f) ? -1 : 1) > 0 || (bad += 1)
        end
    end
    for (labs, flip) in ((0:6, false), ((0, 7:12...), true))
        opp = Dict{Int,Vector{NTuple{2,Int}}}()
        for (t0, l) in zip(S.tri, S.lab), m in 1:3
            l in labs || continue
            # the interface is oriented outwards of the solid: reverse it for the complement
            t = flip && l == 0 ? (t0[1], t0[3], t0[2]) : t0
            push!(get!(opp, t[m], NTuple{2,Int}[]), (t[mod1(m + 1, 3)], t[mod1(m + 2, 3)]))
        end
        for (_, es) in opp                        # the opposite edges of one fan form a single cycle
            nxt = Dict(a => b for (a, b) in es)
            length(nxt) == length(es) || (bad += 1; continue)
            x = es[1][1]; k = 0
            while true
                x = get(nxt, x, 0); k += 1
                (x == es[1][1] || x == 0 || k > length(es)) && break
            end
            (x == es[1][1] && k == length(es)) || (bad += 1)
        end
    end
    return bad
end

# a closed surface without defects: the default construction, then a translated grid and weaker
# snapping (nodes of the grid on the surface give pinched or degenerate facets)
function _tpms_valid_surface(kind, type, level, n, N)
    for snap in (0.25, 0.1, 0.0), shift in (0.0, 0.37 / N)
        S, F = _tpms_surface(kind, type, level, n, N; snap, shift)
        _tp_watertight(S, 0:6) && _tp_watertight(S, (0, 7:12...)) && _tp_defects(S, F) == 0 && return S
    end
    return nothing
end

# range of f over a period (sampled), for the admissible levels
function _tpms_range(kind)
    f = _TPMS_F[kind]
    r = range(0, 2π; length=49)
    vals = [f(x, y, z) for x in r, y in r, z in r]
    return minimum(vals), maximum(vals)
end

function _tpms_level_range(kind, type)
    fmin, fmax = _tpms_range(kind)
    return type === :network ? (fmin, fmax) : (0.0, max(fmax, -fmin))
end

# true if a component of the solid grid nodes wraps around the cell in three independent directions
function _tpms_percolates(kind, type, level, n, N)
    F = _tp_levelfun(kind, type, level, n)
    solid = [F((i / N, j / N, k / N)) <= -1e-9 for i in 0:N-1, j in 0:N-1, k in 0:N-1]
    seen = falses(N, N, N)
    off = Array{NTuple{3,Int}}(undef, N, N, N)
    for s in CartesianIndices(solid)
        (solid[s] && !seen[s]) || continue
        wraps = NTuple{3,Int}[]
        seen[s] = true; off[s] = (0, 0, 0)
        queue = [s]
        while !isempty(queue)
            c = pop!(queue)
            for d in 1:3, sg in (-1, 1)
                q = Tuple(c) .+ ntuple(i -> i == d ? sg : 0, 3)
                w = ntuple(i -> fld(q[i] - 1, N), 3)
                qq = CartesianIndex(ntuple(i -> mod1(q[i], N), 3))
                solid[qq] || continue
                o = off[c] .+ w
                if seen[qq]
                    o == off[qq] || push!(wraps, o .- off[qq])
                else
                    seen[qq] = true; off[qq] = o; push!(queue, qq)
                end
            end
        end
        # rank of the winding vectors (integer, exact in Float64)
        isempty(wraps) && continue
        M = Float64.(reduce(hcat, collect.(unique(wraps))))
        LinearAlgebra.rank(M) == 3 && return true
    end
    return false
end

function _tpms_checked_surface(kind, type, level, n, N, what)
    S = _tpms_valid_surface(kind, type, level, n, N)
    S === nothing && throw(ArgumentError("$what: no valid triangulated surface at this resolution (walls too thin for the " *
                                         "grid, or the surface pinches at the grid nodes); increase resolution or change the density slightly"))
    # mean wall thickness 2V/A of the solid and of the complement, in grid cubes
    cells = _tp_cells(S, N)
    cells >= 2.2 || throw(ArgumentError("$what: the walls are only $(round(cells; digits=2)) grid cubes thick; " *
                                         "increase resolution (cubes per period) or move the density away from 0 and 1"))
    return S
end
_tp_cells(S, N) = (V = _tp_solid_volume(S); 2 * min(V, 1 - V) / _tp_area(S) * N)

# resolution: the given one, or 24 cubes per period raised until the walls are 2.5 cubes thick
function _tpms_resolution(kind, type, level, n, resolution)
    resolution === nothing || return resolution
    res = 24
    S = _tpms_valid_surface(kind, type, level, n, n * res)
    S === nothing && return res
    c = _tp_cells(S, n * res)
    isfinite(c) && c > 0 || return 64
    return c >= 2.5 ? res : min(64, ceil(Int, res * 2.6 / c))
end

function _tpms_check_level(kind, type, level, what)
    lo, hi = _tpms_level_range(kind, type)
    level isa Real && isfinite(level) || throw(ArgumentError("$what: level must be a finite number"))
    lo < level < hi || throw(ArgumentError("$what: the level of a $kind $type must be in ($(round(lo; digits=4)), $(round(hi; digits=4)))"))
    return Float64(level)
end

"""
    tpms_density(kind; type = :sheet, level, periods = 1, resolution = nothing)

Relative density (solid volume fraction) of a triply periodic minimal surface solid in a cubic
cell with `periods` periods per side: `type = :network` is the solid `f ≤ level`, `:sheet` the
wall `|f| ≤ level`, with `f` the level-set function of `kind` (`:gyroid`, `:schwarz_p`,
`:diamond`, `:iwp`, `:neovius`) over one period `[0, 2π]³`. The solid is the faceted surface
of a periodic marching-tetrahedra grid of `resolution` cubes per period (by default 24, raised
until the walls are 2.5 cubes thick), smoothed onto the surface: the one that is meshed. Its
density converges as O(h²). Networks that do not connect across the cell raise an error.
Experimental.
"""
function tpms_density(kind::Symbol; type::Symbol=:sheet, level, periods=1, resolution=nothing)
    n, _ = _tpms_args(kind, type, periods, resolution, "tpms_density")
    c = _tpms_check_level(kind, type, level, "tpms_density")
    res = _tpms_resolution(kind, type, c, n, resolution)
    type === :network && !_tpms_percolates(kind, type, c, n, n * res) &&
        throw(ArgumentError("tpms_density: this $kind network does not connect across the cell (below its pinch-off)"))
    return _tp_solid_volume(_tpms_checked_surface(kind, type, c, n, n * res, "tpms_density"))
end

"""
    tpms_level(kind; type = :sheet, density, periods = 1, resolution = nothing) -> level

Level of [`tpms_density`](@ref) that gives the relative `density`, by bisection on the volume of
the faceted solid (to about 1e-4: the volume jumps by ~1e-5 when the grid topology changes).
Networks below their pinch-off (e.g. Schwarz P below about 0.21) do not connect across the cell,
and walls too thin for the grid cannot be triangulated: both raise an error. Experimental.
"""
tpms_level(kind::Symbol; type::Symbol=:sheet, density, periods=1, resolution=nothing) =
    first(_tpms_level(kind, type, density, periods, resolution, "tpms_level"))

function _tpms_level(kind, type, density, periods, resolution, what)
    n, _ = _tpms_args(kind, type, periods, resolution, what)
    density isa Real && 0 < density < 1 || throw(ArgumentError("$what: density must be in (0, 1)"))
    lo, hi = _tpms_level_range(kind, type)
    function bisect(N)
        a, b = lo + 1e-3, hi - 1e-3
        a0 = a
        c, v = a, NaN
        # start off-centre: level 0 is a degenerate surface for some kinds
        frac = 0.5 + 1 / (2π * 7)
        for _ in 1:60
            c = a + frac * (b - a)
            v = _tp_solid_volume(first(_tpms_surface(kind, type, c, n, N)))
            (abs(v - density) < 1e-4 || b - a < 1e-9 * (hi - lo)) && break
            v < density ? (a = c) : (b = c)
            frac = 0.5
        end
        abs(v - density) < 1e-3 || throw(ArgumentError("$what: density $density is not reachable for a $kind $type " *
                                                       "at this resolution (walls too thin for the grid)"))
        return c
    end
    res = resolution === nothing ? 24 : resolution
    c = try
        bisect(n * res)
    catch e
        # thin walls at the default resolution: try a finer grid
        (e isa ArgumentError && resolution === nothing) || rethrow()
        res = 40
        bisect(n * res)
    end
    if resolution === nothing
        r2 = _tpms_resolution(kind, type, c, n, nothing)
        r2 != res && (res = r2; c = bisect(n * res))
    end
    type === :network && !_tpms_percolates(kind, type, c, n, n * res) &&
        throw(ArgumentError("$what: a $kind network of density $density does not connect across the cell (below its pinch-off); raise the density"))
    # a level whose surface is degenerate: the nearest valid level with the same density (to 1e-3)
    S = _tpms_valid_surface(kind, type, c, n, n * res)
    if S === nothing
        δ = 1e-3 * (hi - lo)
        for k in 1:6, sg in (1, -1)
            c2 = c + sg * k * δ
            S2 = _tpms_valid_surface(kind, type, c2, n, n * res)
            S2 === nothing && continue
            abs(_tp_solid_volume(S2) - density) < 1e-3 && (c = c2; S = S2; break)
        end
    end
    _tpms_checked_surface(kind, type, c, n, n * res, what)
    return c, res
end
