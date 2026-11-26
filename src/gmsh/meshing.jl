# Size fields and mesh generation.

function mesh!(rve::RVE{D}, frame::Frame{D}, part::Partition, opts::MeshOptions) where {D}
    f = gmsh.model.mesh.field
    h = to_frame_length(frame, something(opts.h, minimum(rve.cell.size) / 10))
    fields = Int[]
    bg = f.add("MathEval")
    f.setString(bg, "F", string(h))
    push!(fields, bg)
    listname = D == 3 ? "SurfacesList" : "CurvesList"
    # one Distance + Threshold per distinct (size, distmin, distmax): items with the same
    # characteristic length share a field
    groups = Dict{NTuple{3,Float64},Vector{Int}}()
    for (ip, ph) in enumerate(rve.phases)
        ref = ph.refine
        ref === nothing && continue
        for (ii, item) in enumerate(ph.items)
            pieces = [(D, Int(t)) for (t, own) in part.piece_item if own == (ip, ii)]
            isempty(pieces) && continue
            bnd = gmsh.model.getBoundary(pieces, true, false, false)
            surfs = unique(abs(t) for (d, t) in bnd if d == D - 1)
            isempty(surfs) && continue
            ℓ = to_frame_length(frame, characteristic_length(item))
            hint = ref.h === nothing ? 2π * ℓ / (15 * ref.level) : to_frame_length(frame, ref.h)
            dmin = ref.distmin === nothing ? ℓ / 5 : to_frame_length(frame, ref.distmin)
            dmax = ref.distmax === nothing ? ℓ / 2.5 : to_frame_length(frame, ref.distmax)
            append!(get!(groups, (hint, dmin, dmax), Int[]), surfs)
        end
    end
    for ((hint, dmin, dmax), surfs) in groups
        dist = f.add("Distance")
        f.setNumbers(dist, listname, unique(surfs))
        f.setNumber(dist, "Sampling", 100)
        thr = f.add("Threshold")
        f.setNumber(thr, "InField", dist)
        f.setNumber(thr, "SizeMin", min(hint, h))
        f.setNumber(thr, "SizeMax", h)
        f.setNumber(thr, "DistMin", dmin)
        f.setNumber(thr, "DistMax", max(dmax, dmin * (1 + 1e-9)))
        push!(fields, thr)
    end
    fmin = f.add("Min")
    f.setNumbers(fmin, "FieldsList", fields)
    f.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", opts.curvature)
    try
        gmsh.model.mesh.generate(D)
    catch e
        e isa ErrorException || rethrow()
        error("Gmsh failed to mesh the RVE ($(e.msg)). Common causes: an inclusion crossing or " *
              "touching a face by a tiny amount (sliver), or a global size h much larger than the inclusions.")
    end
    return nothing
end

function mesh_stats(frame::Frame{D}, opts::MeshOptions) where {D}
    ntags, _, _ = gmsh.model.mesh.getNodes()
    _, etags, _ = gmsh.model.mesh.getElements(D)
    tags = isempty(etags) ? Int[] : reduce(vcat, etags)
    isempty(tags) && return MeshStats(length(ntags), 0, opts.order, NaN, NaN, NaN, NaN)
    q = gmsh.model.mesh.getElementQualities(tags, "minSICN")
    emin = gmsh.model.mesh.getElementQualities(tags, "minEdge")
    emax = gmsh.model.mesh.getElementQualities(tags, "maxEdge")
    # called after restore_frame!: node coordinates are physical
    return MeshStats(length(ntags), length(tags), opts.order, minimum(q), sum(q) / length(q),
                     minimum(emin), maximum(emax))
end
