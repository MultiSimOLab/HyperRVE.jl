# Shared test helpers. Tests never assert exact node counts or coordinates.

const OUTDIR = mktempdir()
sphere_vol(r) = 4π * r^3 / 3
unitcell(; kw...) = Cell((1.0, 1.0, 1.0); kw...)
geo(rve) = generate(rve; geometry_only=true)

"Independent O(n^2) check that every node of face_min has a partner on face_max shifted by L."
function brute_periodic_match(fmin, fmax, dir, L; tol=1e-8 * L)
    getn(name) = begin
        t = only(t for (d, t) in gmsh.model.getPhysicalGroups(2) if gmsh.model.getPhysicalName(d, t) == name)
        reshape(gmsh.model.mesh.getNodesForPhysicalGroup(2, t)[2], 3, :)
    end
    pmin, pmax = getn(fmin), getn(fmax)
    size(pmin, 2) > 0 && size(pmin, 2) == size(pmax, 2) || return false
    for i in axes(pmin, 2)
        p = pmin[:, i]
        p[dir] += L
        any(j -> maximum(abs.(pmax[:, j] .- p)) < tol, axes(pmax, 2)) || return false
    end
    return true
end

"Sum of the volumes of the mesh elements of dimension 3."
function mesh_volume()
    _, etags, _ = gmsh.model.mesh.getElements(3)
    return sum(sum(gmsh.model.mesh.getElementQualities(t, "volume")) for t in etags)
end

"Sum of the areas of the mesh elements of dimension 2 (2D cells)."
function mesh_area()
    _, etags, _ = gmsh.model.mesh.getElements(2)
    return sum(sum(gmsh.model.mesh.getElementQualities(t, "volume")) for t in etags)
end

"Independent check of periodic node matching on the sides of a 2D cell."
function brute_periodic_match_2d(fmin, fmax, dir, L; tol=1e-8 * L)
    getn(name) = begin
        t = only(t for (d, t) in gmsh.model.getPhysicalGroups(1) if gmsh.model.getPhysicalName(d, t) == name)
        reshape(gmsh.model.mesh.getNodesForPhysicalGroup(1, t)[2], 3, :)
    end
    pmin, pmax = getn(fmin), getn(fmax)
    size(pmin, 2) > 0 && size(pmin, 2) == size(pmax, 2) || return false
    for i in axes(pmin, 2)
        p = pmin[:, i]
        p[dir] += L
        any(j -> maximum(abs.(pmax[:, j] .- p)) < tol, axes(pmax, 2)) || return false
    end
    return true
end

"""
Connected components of the 2D mesh of physical group `name`, joining elements that share an
edge and identifying the nodes of opposite sides of the periodic cell (points do not connect).
"""
function periodic_components_2d(cell, name; tol=1e-7)
    L, o = cell.size, cell.origin
    tag = only(t for (d, t) in gmsh.model.getPhysicalGroups(2) if gmsh.model.getPhysicalName(d, t) == name)
    key = Dict{Int,NTuple{2,Int}}()
    elems = Vector{Vector{Int}}()
    for e in gmsh.model.getEntitiesForPhysicalGroup(2, tag)
        nt, xyz, _ = gmsh.model.mesh.getNodes(2, e, true)
        for (k, n) in enumerate(nt)
            f = ((xyz[3k-2] - o[1]) / L[1], (xyz[3k-1] - o[2]) / L[2])
            key[Int(n)] = map(x -> round(Int, (x - floor(x + tol)) / tol), f)
        end
        for (t, nn) in zip(gmsh.model.mesh.getElements(2, e)[[1, 3]]...)
            nv = gmsh.model.mesh.getElementProperties(t)[6]       # corner nodes
            np = gmsh.model.mesh.getElementProperties(t)[4]
            for c in 1:np:length(nn)
                push!(elems, [Int(nn[c+k]) for k in 0:nv-1])
            end
        end
    end
    parent = collect(eachindex(elems))
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    owner = Dict{Any,Int}()
    for (i, el) in enumerate(elems), k in eachindex(el)
        a, b = key[el[k]], key[el[mod1(k + 1, length(el))]]
        ed = a < b ? (a, b) : (b, a)
        j = get!(owner, ed, i)
        ri, rj = find(i), find(j)
        ri == rj || (parent[ri] = rj)
    end
    return length(unique(find(i) for i in eachindex(elems)))
end
