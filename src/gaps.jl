# Gap detection and the mesh size callback (see GapRefinement).

# a gap: midpoint, unit axis (across the gap), width, reduced radius R* and number of elements
struct _Gap{D}
    p::NTuple{D,Float64}
    u::NTuple{D,Float64}
    t::Float64
    R::Float64
    face::Int            # 0 between items, d for a face normal to direction d
end

_ball(g::Union{Sphere,Disk}, s) = (g.center .+ s, g.radius)
_ball(g::Translate, s) = _ball(g.child, s .+ g.shift)
_ball(g, s) = nothing

# gaps of width 0 < t < tmax between the balls (with periodic images) and with the cell faces
function _find_gaps(balls::Vector{Tuple{NTuple{D,Float64},Float64}}, lo, ext, per, tmax, faces) where {D}
    rng = ntuple(d -> per[d] ? (-1:1) : (0:0), D)
    found = Dict{NTuple{D,Int},_Gap{D}}()
    key(p) = ntuple(d -> round(Int, (p[d] - lo[d]) / ext[d] * 1e7), D)
    inside(p, w) = all(d -> lo[d] - w < p[d] < lo[d] + ext[d] + w, 1:D)
    add!(g) = (k = key(g.p); (!haskey(found, k) || found[k].t > g.t) && (found[k] = g))
    for (ia, (c0, ra)) in enumerate(balls), ka in Iterators.product(rng...)
        ca = c0 .+ ka .* ext
        if faces
            for d in 1:D, side in (0, 1)
                plane = lo[d] + side * ext[d]
                s = abs(ca[d] - plane) - ra
                0 < s < tmax || continue
                p = ntuple(e -> e == d ? plane + sign(ca[d] - plane) * s / 2 : ca[e], D)
                inside(p, 1e-9) && add!(_Gap{D}(p, ntuple(e -> e == d ? 1.0 : 0.0, D), s, ra, d))
            end
        end
        all(==(0), ka) || continue
        for ib in ia:length(balls), kb in Iterators.product(rng...)
            (ib == ia && all(==(0), kb)) && continue
            cb = balls[ib][1] .+ kb .* ext
            rb = balls[ib][2]
            dv = cb .- ca
            dist = sqrt(sum(abs2, dv))
            t = dist - ra - rb
            0 < t < tmax || continue
            u = dv ./ dist
            p = ca .+ u .* (ra + t / 2)
            R = ra * rb / (ra + rb)
            for k in Iterators.product(rng...)
                q = p .+ k .* ext
                inside(q, sqrt(2R * t)) && add!(_Gap{D}(q, u, t, R, 0))
            end
        end
    end
    return collect(values(found))
end

# ---- size callback: an index of the gaps in cells of the domain, in frame coordinates
struct _GapI
    p::NTuple{3,Float64}
    u::NTuple{3,Float64}
    t::Float64
    twoR::Float64
    n::Float64
    hg::Float64
end
struct _GapIndex
    gaps::Vector{_GapI}
    lo::NTuple{3,Float64}
    ext::NTuple{3,Float64}
    per::NTuple{3,Bool}
    nb::NTuple{3,Int}
    bins::Vector{Vector{Int32}}
    g1::Float64
end
const _GAP_INDEX = Ref{_GapIndex}()

# size of a gap at offset δ from its midpoint: width of the ligament w = t + ρ²/(2R*) over n, growing
# with slope g1 into the inclusions, within a cone from the midpoint
@inline function _gap_law(G::_GapI, δ::NTuple{3,Float64}, g1::Float64)
    d2 = δ[1]^2 + δ[2]^2 + δ[3]^2
    ζ = abs(δ[1] * G.u[1] + δ[2] * G.u[2] + δ[3] * G.u[3])
    w = G.t + max(d2 - ζ^2, 0.0) / G.twoR
    return max(G.hg, min(w / G.n + g1 * max(0.0, ζ - w / 2), G.hg + g1 * sqrt(d2)))
end

function _gap_size(I::_GapIndex, x::NTuple{3,Float64}, lc::Float64)
    b = 0
    stride = 1
    for d in 1:3
        s = (x[d] - I.lo[d]) / I.ext[d]
        I.per[d] && (s -= floor(s))
        k = clamp(floor(Int, s * I.nb[d]), 0, I.nb[d] - 1)
        b += k * stride
        stride *= I.nb[d]
    end
    h = lc
    @inbounds for i in I.bins[b+1]
        G = I.gaps[i]
        δ = ntuple(d -> (v = x[d] - G.p[d]; I.per[d] ? v - I.ext[d] * round(v / I.ext[d]) : v), 3)
        h = min(h, _gap_law(G, δ, I.g1))
    end
    return h
end

# module-level function: its @cfunction is a plain pointer (closures are not supported on every platform)
_gap_callback(dim::Cint, tag::Cint, x::Cdouble, y::Cdouble, z::Cdouble, lc::Cdouble, data::Ptr{Cvoid})::Cdouble =
    _gap_size(_GAP_INDEX[], (x, y, z), lc)

# installs the callback for the gaps of the RVE (frame coordinates); returns the number of gaps
function _install_gaps!(rve::RVE{D}, frame::Frame{D}, opts::MeshOptions, H, hmin) where {D}
    gr = opts.gaps
    lo = to_frame(frame, domain_origin(rve))
    ext = to_frame_length.(Ref(frame), domain_size(rve))
    per = rve.cell.periodic
    balls = Tuple{NTuple{D,Float64},Float64}[]
    unsupported = false
    for k in Iterators.product(ntuple(d -> 0:rve.tiling[d]-1, D)...), ph in rve.phases, g in ph.items
        sh = ntuple(d -> k[d] * rve.cell.size[d], D)
        b = _ball(g, ntuple(_ -> 0.0, D))
        b === nothing && (unsupported = true; continue)
        push!(balls, (Tuple(to_frame(frame, Tuple(b[1] .+ sh))), to_frame_length(frame, b[2])))
    end
    isempty(balls) && (@warn "MeshOptions gaps: no spheres or disks in the RVE; gaps has no effect"; return 0)
    unsupported && @warn "MeshOptions gaps: only the gaps between spheres or disks are refined; other items are ignored"
    perd = Tuple(per[d] && true for d in 1:D)
    gs = _find_gaps(balls, Tuple(lo), Tuple(ext), perd, gr.elements * H, gr.faces)
    g1 = something(gr.growth, opts.growth) - 1
    l3(v, z) = D == 3 ? Tuple(v) : (v[1], v[2], z)
    gaps = _GapI[]
    reach = Float64[]
    for G in gs
        n = G.face > 0 && per[G.face] ? max(1.0, gr.elements / 2) : gr.elements
        hg = max(G.t / n, hmin)
        hg >= H && continue
        push!(gaps, _GapI(l3(G.p, 0.0), l3(G.u, 0.0), G.t, 2G.R, n, hg))
        push!(reach, (H - hg) / g1)
    end
    isempty(gaps) && return 0
    lo3, ext3 = l3(lo, -0.5), l3(ext, 1.0)
    per3 = opts.periodic_sizes ? l3(perd, false) : (false, false, false)
    rmax = maximum(reach)
    nb = ntuple(d -> clamp(floor(Int, 2 * ext3[d] / rmax), 1, 48), 3)
    bins = [Int32[] for _ in 1:prod(nb)]
    for (i, G) in enumerate(gaps)
        r = reach[i]
        rngs = ntuple(3) do d
            a = floor(Int, (G.p[d] - r - lo3[d]) / ext3[d] * nb[d])
            b = floor(Int, (G.p[d] + r - lo3[d]) / ext3[d] * nb[d])
            per3[d] ? unique(mod.(a:min(b, a + nb[d] - 1), nb[d])) : collect(clamp(a, 0, nb[d] - 1):clamp(b, 0, nb[d] - 1))
        end
        for k3 in rngs[3], k2 in rngs[2], k1 in rngs[1]
            push!(bins[1+k1+nb[1]*(k2+nb[2]*k3)], Int32(i))
        end
    end
    _GAP_INDEX[] = _GapIndex(gaps, lo3, ext3, per3, nb, bins, g1)
    cf = @cfunction(_gap_callback, Cdouble, (Cint, Cint, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cvoid}))
    ierr = Ref{Cint}(0)
    ccall((:gmshModelMeshSetSizeCallback, gmsh.lib), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cint}), cf, C_NULL, ierr)
    ierr[] == 0 || error("MeshOptions gaps: Gmsh could not install the size callback")
    return length(gaps)
end

function _remove_gaps!()
    ierr = Ref{Cint}(0)
    ccall((:gmshModelMeshRemoveSizeCallback, gmsh.lib), Cvoid, (Ptr{Cint},), ierr)
    return nothing
end
