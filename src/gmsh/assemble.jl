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
    for (ip, ph) in enumerate(rve.phases), (ii, item) in enumerate(ph.items)
        g = to_frame(frame, item)
        tags = build!(g, ctx)
        blo, bhi = _bbox(tags)
        ax = invariant_axis(g)
        ranges = ntuple(d -> _shift_range(blo[d] - lo[d], bhi[d] - lo[d], ℓ[d], n[d],
                                          rve.cell.periodic[d], ax == d, tol), D)
        shifts = vec(collect(Iterators.product(ranges...)))
        if isempty(shifts)
            @warn "item $ii of phase $(ph.name) lies outside the cell and is ignored"
            occ.remove(tags, true)
            continue
        end
        # note: copies of an item much larger than the cell overlap each other and make the
        # global fragment slow (minutes for a radius of 1.5 cells); clipping or fusing the
        # copies first was tried and is either fragile or slower
        has_zero = false
        for k in shifts
            if all(==(0), k)
                copies = tags
                has_zero = true
            else
                copies = DimTag[(Int32(d), Int32(t)) for (d, t) in occ.copy(tags)]
                occ.translate(copies, _lift3(k .* ℓ)...)
            end
            for c in copies
                push!(tools, c)
                push!(owner, (ip, ii))
            end
        end
        has_zero || occ.remove(tags, true)
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
