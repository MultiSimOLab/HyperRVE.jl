# Size fields and mesh generation.

"""
    _image_shifts(boxes, dmax, lo, ext, periodic)

Lattice shifts `k ≠ 0` (in units of the domain size `ext`) for which some box of `boxes`,
translated by `-k ∘ ext` and grown by `dmax`, meets the domain `[lo, lo + ext]`: for those
shifts, points of the domain lie within `dmax` of a periodic image of the surfaces.
Boxes are `(xmin, ymin, zmin, xmax, ymax, zmax)` tuples; only periodic directions shift.
"""
function _image_shifts(boxes, dmax, lo, ext, periodic::NTuple{D,Bool}) where {D}
    m = ntuple(d -> periodic[d] ? max(1, ceil(Int, dmax / ext[d])) : 0, D)
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
    noeffect = Set{String}()
    for s in ifs, (ip, ii) in s.sides
        ip >= 1 || continue
        ph = rve.phases[ip]
        ref = ph.refine
        hr = ref === nothing ? Inf :
             (ref.h === nothing ? _level_size(S(characteristic_length(ph.items[ii])), ref.level, H) : S(ref.h))
        hi = max(min(hr, ph.h === nothing ? Inf : S(ph.h)), hmin)
        if hi >= H
            ref === nothing || push!(noeffect, ph.name)
            continue
        end
        w = ref === nothing ? 0.0 : (ref.width === nothing ? hi : S(ref.width))
        g = ref === nothing || ref.growth === nothing ? opts.growth : ref.growth
        spec = (hi, w, w + (H - hi) / (g - 1))
        push!(get!(specs, spec, Int[]), Int(s.tag))
    end
    for n in noeffect
        @warn "phase $n: refine has no effect (the interface size is not smaller than the bulk size)"
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
        any(d -> rve.cell.periodic[d] && dmax > ext[d], 1:D) &&
            @warn "the transition of a refined interface is longer than the cell; use a larger growth or a smaller h"
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
    gmsh.option.setNumber("Mesh.MeshSizeMax", opts.hmax === nothing ? H : S(opts.hmax))
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
