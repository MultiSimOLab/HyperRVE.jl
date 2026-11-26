# Periodic images, global fragment and phase assignment.

"""Internal: conformal partition of the domain into phases (entity tags of dimension D)."""
struct Partition
    dim::Int
    names::Vector{String}                       # matrix first (if any), then phases
    tags::Vector{Int}                           # physical tag of each phase
    entities::Vector{Vector{Int32}}
    piece_item::Dict{Int32,Tuple{Int,Int}}      # entity -> (phase index, item index); phase 0 = matrix
end

"""
    _shift_range(a, b, ℓ, n, periodic, invariant, tol)

Lattice shifts `k` (cell lengths `ℓ`) for which the interval `[a, b]` (relative to the
domain origin) overlaps or touches the domain `[0, n ℓ]`. Copies that only touch the domain
are kept on purpose: their footprint is imprinted on the face, so that a body flush with a
periodic face leaves the same imprint on the opposite face. Without periodicity nothing wraps
around and only the `n` tiling copies are considered.
"""
function _shift_range(a, b, ℓ, n, periodic, invariant, tol)
    invariant && return 0:0
    overlaps(k) = (a + k * ℓ < n * ℓ + tol) && (b + k * ℓ > -tol)
    if periodic
        kmin = floor(Int, (-tol - b) / ℓ)
        kmax = ceil(Int, (n * ℓ + tol - a) / ℓ)
        return filter(overlaps, kmin:kmax)
    else
        return filter(overlaps, 0:n-1)
    end
end

function _bbox(tags::Vector{DimTag})
    lo = fill(Inf, 3)
    hi = fill(-Inf, 3)
    for (d, t) in tags
        b = gmsh.model.occ.getBoundingBox(d, t)
        lo .= min.(lo, b[1:3])
        hi .= max.(hi, b[4:6])
    end
    return lo, hi
end

"""
Build `g` (frame coordinates) and its periodic/tiling copies that meet the domain. Returns
`(shift, copies)` pairs; an item outside the domain gives a warning and nothing.
"""
function _item_copies!(g::Geometry{D}, ctx, lo, ℓ, n, periodic, tol, what) where {D}
    occ = gmsh.model.occ
    tags = build!(g, ctx)
    blo, bhi = _bbox(tags)
    ax = invariant_axis(g)
    ranges = ntuple(d -> _shift_range(blo[d] - lo[d], bhi[d] - lo[d], ℓ[d], n[d], periodic[d], ax == d, tol), D)
    shifts = vec(collect(Iterators.product(ranges...)))
    out = Tuple{NTuple{D,Int},Vector{DimTag}}[]
    if isempty(shifts)
        @warn "$what lies outside the cell and is ignored"
        occ.remove(tags, true)
        return out
    end
    # note: copies of an item much larger than the cell overlap each other and make the
    # global fragment slow (minutes for a radius of 1.5 cells)
    has_zero = false
    for k in shifts
        if all(==(0), k)
            push!(out, (k, tags))
            has_zero = true
        else
            copies = DimTag[(Int32(d), Int32(t)) for (d, t) in occ.copy(tags)]
            occ.translate(copies, _lift3(k .* ℓ)...)
            push!(out, (k, copies))
        end
    end
    has_zero || occ.remove(tags, true)
    return out
end

# Identity of a placed primitive, to drop exact duplicates before a fuse (coincident bodies
# are the pathological case of OCC booleans).
_placement_key(::Geometry, shift) = nothing
_rk(x) = round.(x; digits=9) .+ 0.0          # + 0.0 turns -0.0 into 0.0 (isequal in a Set)
_placement_key(g::Sphere, s) = (:sphere, _rk(g.center .+ s), round(g.radius; digits=12))
_placement_key(g::Disk, s) = (:disk, _rk(g.center .+ s), round(g.radius; digits=12))
function _placement_key(g::Union{Cylinder,Capsule}, s)
    p, q = g isa Cylinder ? (g.base, g.base .+ g.axis) : (g.p1, g.p2)
    a, b = _rk(p .+ s), _rk(q .+ s)
    return (nameof(typeof(g)), min(a, b), max(a, b), round(g.radius; digits=12))
end

_mass(tags) = sum(gmsh.model.occ.getMass(d, t) for (d, t) in tags; init=0.0)

function _fuse_or_error(objs::Vector{DimTag}, D, name)
    length(objs) <= 1 && return objs
    vin = _mass(objs)
    body = try
        out, _ = gmsh.model.occ.fuse(objs[1:1], objs[2:end])
        _nonempty(_ents(out, D), "Fuse")
    catch e
        e isa ErrorException || rethrow()
        error("phase $name: OCC could not fuse the items ($(e.msg)); for strut lattices try " *
              "node_radius = 1.05 * radius, or fuse = false")
    end
    # a union cannot be larger than the sum of its parts: OCC returns invalid bodies for some
    # degenerate inputs (e.g. coincident spherical caps of capsules ending at cell corners)
    (_mass(body) <= vin * (1 + 1e-6) && all(t -> gmsh.model.occ.getMass(t...) > 0, body)) ||
        error("phase $name: OCC returned an invalid fused body (volume larger than the sum of the items); " *
              "check for coincident or tangent items, or use fuse = false")
    return body
end

"""
Tools of a phase with `fuse = true`: all items and their images united into one body per
connected component. For a supercell periodic (or not tiled) in every direction, one cell is
fused, clipped to the cell and tiled (the same set, much cheaper than fusing every copy).
"""
function _fused_tools!(rve::RVE{D}, ph::Phase{D}, frame, ctx, lo, ℓ, n, tol) where {D}
    per = rve.cell.periodic
    unit = any(>(1), n) && all(d -> per[d] || n[d] == 1, 1:D)
    if unit
        body = _fused_unit_path!(rve, ph, frame, ctx, lo, ℓ, n, tol)
        body === nothing || return body
    end
    return _fused_copies!(rve, ph, frame, ctx, lo, ℓ, n, tol)
end

"All items and their images for an `m`-tiled domain, duplicates dropped, fused."
function _fused_copies!(rve::RVE{D}, ph::Phase{D}, frame, ctx, lo, ℓ, m, tol) where {D}
    occ = gmsh.model.occ
    copies = DimTag[]
    seen = Set{Any}()
    for (ii, item) in enumerate(ph.items)
        g = to_frame(frame, item)
        for (k, c) in _item_copies!(g, ctx, lo, ℓ, m, rve.cell.periodic, tol, "item $ii of phase $(ph.name)")
            key = _placement_key(g, k .* ℓ)
            if key !== nothing && key in seen
                occ.remove(c, true)
            else
                append!(copies, c)
                key === nothing || push!(seen, key)
            end
        end
    end
    isempty(copies) && return copies
    return _fuse_or_error(copies, D, ph.name)
end

"""
Fuse one cell, clip it to the cell and tile it. Returns `nothing` (after removing what it
built) if OCC cannot fuse the tiles, so that the caller falls back to fusing every copy.
"""
function _fused_unit_path!(rve::RVE{D}, ph::Phase{D}, frame, ctx, lo, ℓ, n, tol) where {D}
    occ = gmsh.model.occ
    body = _fused_copies!(rve, ph, frame, ctx, lo, ℓ, ntuple(_ -> 1, D), tol)
    isempty(body) && return body
    cellbox = D == 3 ? occ.addBox(lo..., ℓ...) : occ.addRectangle(lo[1], lo[2], 0.0, ℓ...)
    body = try
        _ents(occ.intersect(body, [(D, cellbox)])[1], D)
    catch e
        e isa ErrorException || rethrow()
        error("phase $(ph.name): OCC could not clip the fused cell ($(e.msg)); use fuse = false")
    end
    isempty(body) && return body
    _mass(body) <= prod(ℓ) * (1 + 1e-6) && all(t -> occ.getMass(t...) > 0, body) ||
        error("phase $(ph.name): OCC returned an invalid fused cell (volume larger than the cell); " *
              "check for coincident or tangent items, or use fuse = false")
    tiles = copy(body)
    for k in Iterators.product(ntuple(d -> 0:n[d]-1, D)...)
        all(==(0), k) && continue
        c = DimTag[(Int32(d), Int32(t)) for (d, t) in occ.copy(body)]
        occ.translate(c, _lift3(k .* ℓ)...)
        append!(tiles, c)
    end
    vin = _mass(tiles)
    try
        out, _ = occ.fuse(tiles[1:1], tiles[2:end])
        fused = _ents(out, D)
        _mass(fused) <= vin * (1 + 1e-6) && return fused
        occ.remove(fused, true)
    catch e
        e isa ErrorException || rethrow()
        occ.remove(filter(t -> t in Set(occ.getEntities(D)), tiles), true)
    end
    return nothing
end

function assemble!(rve::RVE{D}, frame::Frame{D}) where {D}
    occ = gmsh.model.occ
    tol = geomtol(frame)
    ℓ = to_frame_length.(Ref(frame), rve.cell.size)
    n = rve.tiling
    lo = to_frame(frame, domain_origin(rve))
    ext = n .* ℓ
    ctx = BuildContext(_lift3(lo .+ ext ./ 2), _norm(ext))

    tools = DimTag[]
    owner = Tuple{Int,Int}[]
    for (ip, ph) in enumerate(rve.phases)
        if ph.fuse
            body = _fused_tools!(rve, ph, frame, ctx, lo, ℓ, n, tol)
            istar = argmin([characteristic_length(it) for it in ph.items])
            append!(tools, body)
            append!(owner, fill((ip, istar), length(body)))
        else
            for (ii, item) in enumerate(ph.items)
                for (_, copies) in _item_copies!(to_frame(frame, item), ctx, lo, ℓ, n, rve.cell.periodic, tol,
                                                 "item $ii of phase $(ph.name)")
                    append!(tools, copies)
                    append!(owner, fill((ip, ii), length(copies)))
                end
            end
        end
    end

    box = D == 3 ? occ.addBox(lo..., ext...) : occ.addRectangle(lo[1], lo[2], 0.0, ext...)
    piece_phase = Dict{Int32,Int}()
    piece_item = Dict{Int32,Tuple{Int,Int}}()
    inside = Set{DimTag}()
    if isempty(tools)
        push!(inside, (Int32(D), Int32(box)))
        out = DimTag[(Int32(D), Int32(box))]
    else
        out, maps = occ.fragment([(D, box)], tools)
        out = DimTag[(Int32(d), Int32(t)) for (d, t) in out]
        for c in maps[1]
            push!(inside, (Int32(c[1]), Int32(c[2])))
        end
        for (j, children) in enumerate(maps[2:end])
            ip, ii = owner[j]
            for c in children
                key = (Int32(c[1]), Int32(c[2]))
                key in inside || continue
                prev = get(piece_phase, key[2], 0)
                if prev != 0 && prev != ip && rve.overlaps == :error
                    v = from_frame_measure(frame, occ.getMass(c[1], c[2]))
                    a, b = rve.phases[min(prev, ip)].name, rve.phases[max(prev, ip)].name
                    error("phases $a and $b overlap ($(D == 3 ? "volume" : "area") ≈ $(round(v; sigdigits=3))); " *
                          "use overlaps = :paint or separate them")
                end
                if ip > prev
                    piece_phase[key[2]] = ip
                    piece_item[key[2]] = (ip, ii)
                end
            end
        end
    end
    outside = [e for e in out if e[1] == D && !(e in inside)]
    isempty(outside) || occ.remove(outside, true)

    matrix_pieces = Int32[t for (d, t) in inside if get(piece_phase, t, 0) == 0]
    if rve.matrix === nothing && !isempty(matrix_pieces)
        occ.remove([(Int32(D), t) for t in matrix_pieces], true)
    end
    occ.synchronize()

    names = String[]
    ptags = Int[]
    ents = Vector{Int32}[]
    if rve.matrix !== nothing
        push!(names, rve.matrix); push!(ptags, 1); push!(ents, sort(matrix_pieces))
        for t in matrix_pieces
            piece_item[t] = (0, 0)
        end
    end
    for (ip, ph) in enumerate(rve.phases)
        push!(names, ph.name)
        push!(ptags, 1 + ip)
        push!(ents, sort([t for (t, p) in piece_phase if p == ip]))
    end
    for (name, e) in zip(names, ents)
        isempty(e) || continue
        if name == rve.matrix
            @warn "the matrix $name has no volume (the phases fill the cell); consider matrix = nothing"
        else
            @warn "phase $name has no volume inside the cell (covered by later phases or outside)"
        end
    end
    return Partition(D, names, ptags, ents, piece_item)
end
