# Periodic constraints on opposite faces (curves in 2D) and their verification after meshing.

function _translation(d, L)
    T = [1.0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    T[4*d] = L
    return T
end

function set_periodicity!(rve::RVE{D}, frame::Frame{D}, on=_boundary_entities(rve, frame)) where {D}
    any(rve.cell.periodic) || return nothing
    # pairing compares exact (padded) bounding boxes of an entity and its translate: the
    # padding is the same on both sides, but can reach ~1e-5 for B-spline cuts
    tol = 1000 * geomtol(frame)
    ext = to_frame_length.(Ref(frame), domain_size(rve))
    for d in 1:D
        rve.cell.periodic[d] || continue
        masters = sort([t for ((dim, t), ps) in on.planes if dim == D - 1 && (d, 1) in ps])
        slaves = sort([t for ((dim, t), ps) in on.planes if dim == D - 1 && (d, 2) in ps])
        length(masters) == length(slaves) ||
            throw(PeriodicityError(d, "the $(length(masters)) entities on $(_face_name(d, 1)) and the " *
                                      "$(length(slaves)) on $(_face_name(d, 2)) cannot be paired; check for " *
                                      "inclusions tangent or nearly tangent to a face"))
        used = Set{Int32}()
        for m in masters
            bm = collect(on.bbox[(Int32(D - 1), m)])
            bm[d] += ext[d]
            bm[3+d] += ext[d]
            cands = [s for s in slaves if !(s in used) &&
                     maximum(abs.(collect(on.bbox[(Int32(D - 1), s)]) .- bm)) < tol]
            if length(cands) > 1
                cm = collect(gmsh.model.occ.getCenterOfMass(D - 1, m))
                cm[d] += ext[d]
                dist(s) = maximum(abs.(collect(gmsh.model.occ.getCenterOfMass(D - 1, s)) .- cm))
                cands = [cands[argmin(dist.(cands))]]
            end
            isempty(cands) &&
                throw(PeriodicityError(d, "entity $m on $(_face_name(d, 1)) has no periodic partner on " *
                                          "$(_face_name(d, 2)); check for inclusions tangent or nearly " *
                                          "tangent to a face"))
            s = only(cands)
            push!(used, s)
            try
                gmsh.model.mesh.setPeriodic(D - 1, [s], [m], _translation(d, ext[d]))
            catch e
                e isa ErrorException || rethrow()
                throw(PeriodicityError(d, "Gmsh could not match the boundaries of entity $m on " *
                                          "$(_face_name(d, 1)) and entity $s on $(_face_name(d, 2)) " *
                                          "($(sprint(showerror, e))); this usually means an inclusion " *
                                          "tangent to a face from inside the cell"))
            end
        end
    end
    return nothing
end

"""
    _physical_periodic_affine!(path, frame)

The periodic links of a `.msh` file store the affine map between master and slave entities,
and Gmsh applies it to the slave nodes when the file is read. The model is meshed in the
normalised frame, so the stored maps are in frame units: rewrite them in physical units,
`A_phys = T A T⁻¹` with `T(x) = origin + scale x`.
"""
function _physical_periodic_affine!(path::AbstractString, frame::Frame)
    _isidentity(frame) && return nothing
    any(e -> endswith(lowercase(path), e), (".msh", ".msh4")) || return nothing
    lines = readlines(path)
    i = findfirst(==("\$Periodic"), lines)
    i === nothing && return nothing
    s = frame.scale
    o = collect(_lift3(frame.origin))
    T = [s 0 0 o[1]; 0 s 0 o[2]; 0 0 s o[3]; 0 0 0 1]
    Ti = [1/s 0 0 -o[1]/s; 0 1/s 0 -o[2]/s; 0 0 1/s -o[3]/s; 0 0 0 1]
    n = parse(Int, lines[i+1])
    j = i + 2
    for _ in 1:n
        j += 1                                   # entityDim entityTag entityTagMaster
        f = split(lines[j])
        na = parse(Int, f[1])
        if na == 16
            A = permutedims(reshape(parse.(Float64, f[2:17]), 4, 4))
            B = T * A * Ti
            lines[j] = "16 " * join((repr(x) for x in vec(permutedims(B))), " ")
        elseif na != 0
            @warn "unexpected periodic affine map with $na values in $path; left unchanged"
        end
        j += 1
        np = parse(Int, lines[j])
        j += np + 1
    end
    open(path, "w") do io
        for l in lines
            println(io, l)
        end
    end
    return nothing
end

function _group_nodes(dim, name)
    for (d, t) in gmsh.model.getPhysicalGroups(dim)
        if gmsh.model.getPhysicalName(d, t) == name
            _, c = gmsh.model.mesh.getNodesForPhysicalGroup(d, t)
            return reshape(c, 3, :)
        end
    end
    return zeros(3, 0)
end

function verify_periodicity(rve::RVE{D}, frame::Frame{D}) where {D}
    checks = PeriodicityCheck[]
    L = domain_size(rve)
    lref = maximum(L)
    tol = PERIODIC_TOL * lref
    q = 1e-6 * lref
    for d in 1:D
        rve.cell.periodic[d] || continue
        pmin = _group_nodes(D - 1, _face_name(d, 1))
        pmax = _group_nodes(D - 1, _face_name(d, 2))
        others = [i for i in 1:3 if i != d]
        key(p) = Tuple(round.(Int, p[others] ./ q))
        bins = Dict{NTuple{2,Int},Vector{Int}}()
        for j in axes(pmax, 2)
            push!(get!(bins, key(pmax[:, j]), Int[]), j)
        end
        maxerr = 0.0
        ok = size(pmin, 2) == size(pmax, 2)
        for i in axes(pmin, 2)
            p = pmin[:, i]
            p[d] += L[d]
            k = key(p)
            best = Inf
            for a in -1:1, b in -1:1, j in get(bins, (k[1] + a, k[2] + b), Int[])
                best = min(best, maximum(abs.(pmax[:, j] .- p)))
            end
            maxerr = max(maxerr, best)
            best <= tol || (ok = false)
        end
        # with a matrix the cell faces always exist: finding no nodes means they were not
        # recognised, which must not pass silently
        rve.matrix !== nothing && size(pmin, 2) == 0 && (ok = false)
        push!(checks, PeriodicityCheck(d, size(pmin, 2), maxerr / lref, ok))
        ok || throw(PeriodicityError(d, "the meshes of $(_face_name(d, 1)) and $(_face_name(d, 2)) do not " *
                                        "match ($(size(pmin, 2)) vs $(size(pmax, 2)) nodes, max mismatch " *
                                        "$(maxerr / lref) relative)"))
    end
    return checks
end
