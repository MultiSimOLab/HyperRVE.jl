# Physical groups: phases (dimension D), cell faces/edges/corners (lower dimensions).

const _SIDE = ("min", "max")

_face_name(d, s) = "xyz"[d] * _SIDE[s]
_face_tag(d, s) = 2 * (d - 1) + s

"Planes (d, side) on which an entity lies, from its bounding box."
function _planes(dim, tag, lo, hi, D, tol)
    b = gmsh.model.getBoundingBox(dim, tag)
    ps = Tuple{Int,Int}[]
    for d in 1:D
        (abs(b[d] - lo[d]) < tol && abs(b[3+d] - lo[d]) < tol) && push!(ps, (d, 1))
        (abs(b[d] - hi[d]) < tol && abs(b[3+d] - hi[d]) < tol) && push!(ps, (d, 2))
    end
    return ps
end

"""Lower-dimensional entities on the cell boundary: Dict (dim, tag) => planes."""
function _boundary_entities(rve::RVE{D}, frame::Frame{D}) where {D}
    lo = to_frame(frame, domain_origin(rve))
    hi = lo .+ to_frame_length.(Ref(frame), domain_size(rve))
    tol = geomtol(frame)
    on = Dict{DimTag,Vector{Tuple{Int,Int}}}()
    for dim in 0:D-1, (d, t) in gmsh.model.getEntities(dim)
        ps = _planes(d, t, lo, hi, D, tol)
        isempty(ps) || (on[(Int32(d), Int32(t))] = ps)
    end
    return on
end

function _add_group!(groups, dim, ents, tag, name)
    isempty(ents) && return
    e = sort!(unique(Int.(ents)))
    gmsh.model.addPhysicalGroup(dim, e, tag, name)
    push!(groups, GroupInfo(name, dim, tag, e))
end

function add_physical_groups!(rve::RVE{D}, frame::Frame{D}, part::Partition) where {D}
    groups = GroupInfo[]
    for (name, tag, e) in zip(part.names, part.tags, part.entities)
        _add_group!(groups, D, e, tag, name)
    end
    on = _boundary_entities(rve, frame)
    sel(dim, pred) = [t for ((d, t), ps) in on if d == dim && pred(ps)]
    # closed faces, in every lower dimension
    for d in 1:D, s in 1:2, dim in D-1:-1:0
        _add_group!(groups, dim, sel(dim, ps -> (d, s) in ps), _face_tag(d, s), _face_name(d, s))
    end
    # all corner points
    _add_group!(groups, 0, sel(0, ps -> length(unique(first.(ps))) == D), 7, "corners")
    # whole boundary
    for dim in D-1:-1:0
        _add_group!(groups, dim, sel(dim, _ -> true), 8, "boundary")
    end
    # edges (3D): parallel to x, then y, then z
    if D == 3
        tag = 11
        for a in 1:3
            d1, d2 = Tuple(setdiff(1:3, a))
            for s1 in 1:2, s2 in 1:2
                name = "edge_" * _face_name(d1, s1) * "_" * _face_name(d2, s2)
                for dim in 1:-1:0
                    _add_group!(groups, dim, sel(dim, ps -> (d1, s1) in ps && (d2, s2) in ps), tag, name)
                end
                tag += 1
            end
        end
    end
    # individual corners: x fastest, z slowest
    tag = 31
    for idx in Iterators.product(ntuple(_ -> 1:2, D)...)
        name = "corner_" * join((_face_name(d, idx[d]) for d in 1:D), "_")
        _add_group!(groups, 0, sel(0, ps -> all(d -> (d, idx[d]) in ps, 1:D)), tag, name)
        tag += 1
    end
    return groups
end
