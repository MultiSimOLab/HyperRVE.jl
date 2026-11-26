# Physical groups: phases (dimension D), cell faces/edges/corners (lower dimensions).

const _SIDE = ("min", "max")

_face_name(d, s) = "xyz"[d] * _SIDE[s]
_face_tag(d, s) = 2 * (d - 1) + s

"Planes (d, side) on which a bounding box lies (both extremes within `tol` of the plane)."
function _planes(b, lo, hi, D, tol)
    ps = Tuple{Int,Int}[]
    for d in 1:D
        (abs(b[d] - lo[d]) < tol && abs(b[3+d] - lo[d]) < tol) && push!(ps, (d, 1))
        (abs(b[d] - hi[d]) < tol && abs(b[3+d] - hi[d]) < tol) && push!(ps, (d, 2))
    end
    return ps
end

"""Lower-dimensional entities on the cell boundary, with their planes and exact bounding boxes."""
struct BoundaryEntities
    planes::Dict{DimTag,Vector{Tuple{Int,Int}}}
    bbox::Dict{DimTag,NTuple{6,Float64}}
end

# Exact OCC bounding boxes are padded: ~1e-7 for analytic entities, but up to ~1e-5 for planar
# faces and curves cut by B-spline surfaces (dilated spheres = ellipsoids, non-uniform
# scalings). Bounding boxes of the tessellation are exact in the normal direction of a planar
# entity, but tessellating every entity is slow. So: select candidates with the exact boxes and
# a loose tolerance, then confirm them with the tessellated boxes and the strict tolerance.
const _CANDIDATE_TOL = 1e-3

function _boundary_entities(rve::RVE{D}, frame::Frame{D}) where {D}
    lo = to_frame(frame, domain_origin(rve))
    hi = lo .+ to_frame_length.(Ref(frame), domain_size(rve))
    tol = geomtol(frame)
    loose = _CANDIDATE_TOL * maximum(hi .- lo)
    cands = DimTag[]
    exact = Dict{DimTag,NTuple{6,Float64}}()
    for dim in 0:D-1, (d, t) in gmsh.model.getEntities(dim)
        b = Tuple(gmsh.model.getBoundingBox(d, t))
        if !isempty(_planes(b, lo, hi, D, loose))
            push!(cands, (Int32(d), Int32(t)))
            exact[(Int32(d), Int32(t))] = b
        end
    end
    planes = Dict{DimTag,Vector{Tuple{Int,Int}}}()
    bbox = Dict{DimTag,NTuple{6,Float64}}()
    stl = gmsh.option.getNumber("Geometry.OCCBoundsUseStl")
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
    try
        for e in cands
            ps = _planes(gmsh.model.getBoundingBox(e...), lo, hi, D, tol)
            isempty(ps) && continue
            planes[e] = ps
            bbox[e] = exact[e]     # exact box: same padding for an entity and its translate
        end
    finally
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", stl)
    end
    return BoundaryEntities(planes, bbox)
end

function _add_group!(groups, dim, ents, tag, name)
    isempty(ents) && return
    e = sort!(unique(Int.(ents)))
    gmsh.model.addPhysicalGroup(dim, e, tag, name)
    push!(groups, GroupInfo(name, dim, tag, e))
end

function add_physical_groups!(rve::RVE{D}, frame::Frame{D}, part::Partition,
                              on=_boundary_entities(rve, frame)) where {D}
    groups = GroupInfo[]
    for (name, tag, e) in zip(part.names, part.tags, part.entities)
        _add_group!(groups, D, e, tag, name)
    end
    sel(dim, pred) = [t for ((d, t), ps) in on.planes if d == dim && pred(ps)]
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
