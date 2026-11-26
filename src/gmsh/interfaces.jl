# True interfaces: (D-1)-entities not on a cell face whose adjacent volumes belong to
# different phases (or with a single adjacent volume: free surface, side "void").

struct InterfaceSurface
    tag::Int32
    sides::NTuple{2,Tuple{Int,Int}}   # (phase index, item index); phase 0 = matrix, (-1, 0) = void
end

function classify_interfaces(part::Partition, on::BoundaryEntities)
    D = part.dim
    out = InterfaceSurface[]
    for (d, t) in gmsh.model.getEntities(D - 1)
        haskey(on.planes, (Int32(d), Int32(t))) && continue
        up, _ = gmsh.model.getAdjacencies(D - 1, t)
        own = [part.piece_item[Int32(v)] for v in up if haskey(part.piece_item, Int32(v))]
        if length(own) == 2
            own[1][1] != own[2][1] && push!(out, InterfaceSurface(Int32(t), (own[1], own[2])))
        elseif length(own) == 1
            push!(out, InterfaceSurface(Int32(t), (own[1], (-1, 0))))
        end
    end
    return out
end


"""
Interface physical groups `interface_<a>_<b>` (phases in partition order, `void` last) in
dimensions D-1 … 0 (closed sets, like the face groups), plus their union `interfaces`.
Tags depend only on the phase list: union 99, pair (i, j) -> 100 + (i - 1)(N + 1) + j.
"""
function interface_groups!(groups, rve::RVE{D}, frame::Frame{D}, part::Partition,
                           ifs::Vector{InterfaceSurface}; create::Bool=true) where {D}
    names = vcat(part.names, "void")
    N = length(part.names)
    offset = rve.matrix === nothing ? 0 : 1
    slot(p) = p == -1 ? N + 1 : p + offset
    pairs = Dict{Tuple{Int,Int},Vector{Int}}()
    for s in ifs
        i, j = slot(s.sides[1][1]), slot(s.sides[2][1])
        i, j = minmax(i, j)
        push!(get!(pairs, (i, j), Int[]), Int(s.tag))
    end
    infos = InterfaceInfo[]
    # closed sets: the surfaces and their bounding curves and points
    closure(surfs, dim) = dim == D - 1 ? surfs :
        unique(abs(t) for (d, t) in gmsh.model.getBoundary([(D - 1, s) for s in surfs], false, false, dim == 0)
               if d == dim)
    allsurfs = Int[]
    for (i, j) in sort!(collect(keys(pairs)))
        surfs = sort!(unique(pairs[(i, j)]))
        append!(allsurfs, surfs)
        name = "interface_" * names[i] * "_" * names[j]
        tag = 100 + (i - 1) * (N + 1) + j
        area = sum(from_frame_measure(frame, gmsh.model.occ.getMass(D - 1, s), D - 1) for s in surfs)
        if create
            for dim in D-1:-1:0
                _add_group!(groups, dim, closure(surfs, dim), tag, name)
            end
        end
        push!(infos, InterfaceInfo(name, (names[i], names[j]), create ? tag : 0, surfs, area))
    end
    if create && !isempty(allsurfs)
        u = sort!(unique(allsurfs))
        for dim in D-1:-1:0
            _add_group!(groups, dim, closure(u, dim), 99, "interfaces")
        end
    end
    return infos
end
