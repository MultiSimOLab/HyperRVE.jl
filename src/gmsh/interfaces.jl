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
