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
