# Size fields and mesh generation.

"""
    _image_shifts(boxes, dmax, lo, ext, periodic)

Lattice shifts `k ≠ 0` (in units of the domain size `ext`) for which some box of `boxes`,
translated by `-k ∘ ext` and grown by `dmax`, meets the domain `[lo, lo + ext]`: for those
shifts, points of the domain lie within `dmax` of a periodic image of the surfaces.
Boxes are `(xmin, ymin, zmin, xmax, ymax, zmax)` tuples; only periodic directions shift.
"""
function _image_shifts(boxes, dmax, lo, ext, periodic::NTuple{D,Bool}) where {D}
    # Σ lies inside the domain, so for a point of the domain the image k + sign(k) e_d is never
    # closer than k: shifts beyond the nearest neighbours are not needed, whatever dmax is
    m = ntuple(d -> periodic[d] ? 1 : 0, D)
    ks = Set{NTuple{D,Int}}()
    for k in Iterators.product(ntuple(d -> -m[d]:m[d], D)...)
        all(==(0), k) && continue
        for b in boxes
            if all(d -> b[d] - k[d] * ext[d] - dmax < lo[d] + ext[d] && b[3+d] - k[d] * ext[d] + dmax > lo[d], 1:D)
                push!(ks, k)
                break
            end
        end
    end
    return sort!(collect(ks))
end

# Size of a refinement given by `level` (2024 rule), rounded down to the ladder H·2^(-j/4)
# so that polydisperse items share fields.
function _level_size(ℓ, level, H)
    hl = 2π * ℓ / (15 * level)
    hl >= H && return hl
    j = ceil(Int, -4 * log2(hl / H))
    return H * 2.0^(-j / 4)
end

function mesh!(rve::RVE{D}, frame::Frame{D}, part::Partition, ifs::Vector{InterfaceSurface},
               opts::MeshOptions) where {D}
    f = gmsh.model.mesh.field
    S(x) = to_frame_length(frame, x)
    h = S(something(opts.h, minimum(rve.cell.size) / 10))
    # bulk sizes per phase (matrix: global h)
    bulk(ip) = ip >= 1 && rve.phases[ip].h !== nothing ? S(rve.phases[ip].h) : h
    sizes = Dict{Float64,Vector{Int}}()
    for (t, (ip, _)) in part.piece_item
        push!(get!(sizes, bulk(ip), Int[]), Int(t))
    end
    H = isempty(sizes) ? h : maximum(keys(sizes))
    fields = Int[]
    bg = f.add("MathEval")
    f.setString(bg, "F", string(H))
    push!(fields, bg)
    if length(sizes) > 1
        for (b, vols) in sizes
            c = f.add("Constant")
            f.setNumber(c, "VIn", b)
            f.setNumbers(c, D == 3 ? "VolumesList" : "SurfacesList", vols)
            f.setNumber(c, "IncludeBoundary", 1)
            push!(fields, c)
        end
    end
    # graded cones around the interfaces of refined phases (and of phases finer than H)
    hmin = opts.hmin === nothing ? 0.0 : S(opts.hmin)
    specs = Dict{NTuple{3,Float64},Vector{Int}}()
    effective = Set{String}()
    refined = Set{String}()
    for s in ifs, (side, (ip, ii)) in enumerate(s.sides)
        ip >= 1 || continue
        ph = rve.phases[ip]
        ref = ph.refine
        hr = ref === nothing ? Inf :
             (ref.h === nothing ? _level_size(S(characteristic_length(ph.items[ii])), ref.level, H) : S(ref.h))
        if ref !== nothing
            push!(refined, ph.name)
            other = s.sides[3-side][1]
            # the refinement matters if it is finer than the bulk size on either side
            hr < max(bulk(ip), other >= 0 ? bulk(other) : 0.0) && push!(effective, ph.name)
        end
        hi = max(min(hr, ph.h === nothing ? Inf : S(ph.h)), hmin)
        hi >= H && continue
        w = ref === nothing ? 0.0 : (ref.width === nothing ? hi : S(ref.width))
        g = ref === nothing || ref.growth === nothing ? opts.growth : ref.growth
        spec = (hi, w, w + (H - hi) / (g - 1))
        push!(get!(specs, spec, Int[]), Int(s.tag))
    end
    for n in setdiff(refined, effective)
        @warn "phase $n: refine has no effect (the interface size is not smaller than the bulk sizes on either side)"
    end
    listname = D == 3 ? "SurfacesList" : "CurvesList"
    lo = to_frame(frame, domain_origin(rve))
    ext = to_frame_length.(Ref(frame), domain_size(rve))
    for ((hi, w, dmax), surfs) in specs
        surfs = unique(surfs)
        boxes = [Tuple(gmsh.model.getBoundingBox(D - 1, t)) for t in surfs]
        diag = maximum(b -> sqrt(sum(abs2, (b[4] - b[1], b[5] - b[2], b[6] - b[3]))), boxes)
        dist = f.add("Distance")
        f.setNumbers(dist, listname, surfs)
        f.setNumber(dist, "Sampling", clamp(ceil(Int, 2 * diag / hi), 20, 300))
        thr = f.add("Threshold")
        f.setNumber(thr, "InField", dist)
        f.setNumber(thr, "SizeMin", hi)
        f.setNumber(thr, "SizeMax", H)
        f.setNumber(thr, "DistMin", w)
        f.setNumber(thr, "DistMax", max(dmax, w * (1 + 1e-9) + 1e-12))
        push!(fields, thr)
        opts.periodic_sizes || continue
        # periodic size field: the same field evaluated at lattice-shifted points
        shifts = _image_shifts(boxes, dmax, lo, ext, rve.cell.periodic)
        for k in shifts
            p = f.add("Param")
            f.setNumber(p, "InField", thr)
            sh = _lift3(k .* ext)
            f.setString(p, "FX", @sprintf("x+(%.17f)", sh[1]))
            f.setString(p, "FY", @sprintf("y+(%.17f)", sh[2]))
            f.setString(p, "FZ", @sprintf("z+(%.17f)", sh[3]))
            push!(fields, p)
        end
    end
    fmin = f.add("Min")
    f.setNumbers(fmin, "FieldsList", fields)
    f.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", opts.curvature)
    gmsh.option.setNumber("Mesh.MeshSizeMin", hmin)
    # hmin is a hard lower bound: it wins over the bulk sizes as well
    gmsh.option.setNumber("Mesh.MeshSizeMax", max(opts.hmax === nothing ? H : S(opts.hmax), hmin))
    gmsh.option.setNumber("Mesh.Algorithm", _ALG2D[opts.surface_algorithm])
    gmsh.option.setNumber("Mesh.Algorithm3D", _ALG3D[opts.algorithm])
    gmsh.option.setNumber("Mesh.Optimize", opts.optimize == :none ? 0 : 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", opts.optimize == :netgen ? 1 : 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)          # generation is always linear
    try
        gmsh.model.mesh.generate(D)
        if opts.order == 2
            # before restore_frame!: curved mid-nodes are projected on the (normalised) CAD
            gmsh.option.setNumber("Mesh.SecondOrderLinear", opts.curved ? 0 : 1)
            gmsh.option.setNumber("Mesh.HighOrderOptimize", 0)
            gmsh.model.mesh.setOrder(2)
            opts.optimize_curved == :none || gmsh.model.mesh.optimize(_HO[opts.optimize_curved])
        end
    catch e
        e isa ErrorException || rethrow()
        alt = opts.algorithm == :delaunay ? "" : ", or try algorithm = :delaunay"
        error("Gmsh failed to mesh the RVE ($(e.msg)). Common causes: an inclusion crossing or " *
              "touching a face by a tiny amount (sliver), or a global size h much larger than the inclusions$alt.")
    end
    return nothing
end

function mesh_stats(frame::Frame{D}, opts::MeshOptions, part::Partition) where {D}
    ntags, _, _ = gmsh.model.mesh.getNodes()
    nan3 = (NaN, NaN, NaN)
    # elements per phase (one getElements per entity of dimension D)
    tags = Int[]
    owner = Int[]                       # index into part.names
    for (k, ents) in enumerate(part.entities), e in ents
        _, et, _ = gmsh.model.mesh.getElements(D, e)
        for t in et
            append!(tags, t)
            append!(owner, fill(k, length(t)))
        end
    end
    isempty(tags) && return MeshStats(length(ntags), 0, opts.order, NaN, NaN, NaN, NaN, ntuple(_ -> 0, 10), 0,
                                      NaN, nan3, "", Pair{String,Float64}[])
    q = gmsh.model.mesh.getElementQualities(tags, "minSICN")
    emin = gmsh.model.mesh.getElementQualities(tags, "minEdge")
    emax = gmsh.model.mesh.getElementQualities(tags, "maxEdge")
    sj = opts.order == 2 ? gmsh.model.mesh.getElementQualities(tags, "minSJ") : nothing
    hist = zeros(Int, 10)
    for x in q
        hist[clamp(floor(Int, 10x) + 1, 1, 10)] += 1
    end
    ninv = count(<=(0), sj === nothing ? q : sj)
    w = argmin(q)
    _, nodes, _, _ = gmsh.model.mesh.getElement(tags[w])
    pts = [gmsh.model.mesh.getNode(n)[1] for n in nodes]
    wp = Tuple(sum(pts) ./ length(pts))
    pq = Pair{String,Float64}[]
    for (k, name) in enumerate(part.names)
        sel = q[owner.==k]
        isempty(sel) || push!(pq, name => minimum(sel))
    end
    # called after restore_frame!: node coordinates are physical
    return MeshStats(length(ntags), length(tags), opts.order, minimum(q), sum(q) / length(q),
                     minimum(emin), maximum(emax), Tuple(hist), ninv,
                     sj === nothing ? NaN : minimum(sj), (wp[1], wp[2], wp[3]), part.names[owner[w]], pq)
end
