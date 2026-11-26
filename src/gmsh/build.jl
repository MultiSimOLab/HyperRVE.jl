# Geometry tree -> OCC entities (no synchronize). Returns the D-dimensional entities.

struct BuildContext
    center::Vec{3}     # centre of the meshed domain (frame coordinates, z = 0 in 2D)
    diag::Float64      # length of the domain diagonal
end

_ents(out, D) = DimTag[(Int32(d), Int32(t)) for (d, t) in out if d == D]

function _nonempty(tags, what)
    isempty(tags) && error("boolean operation $what produced an empty geometry")
    return tags
end

build!(g::Sphere, ::BuildContext) = DimTag[(3, gmsh.model.occ.addSphere(g.center..., g.radius))]

function build!(g::Ellipsoid, ::BuildContext)
    t = DimTag[(3, gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0))]
    gmsh.model.occ.dilate(t, 0.0, 0.0, 0.0, g.radii...)
    gmsh.model.occ.translate(t, g.center...)
    return t
end

build!(g::Box, ::BuildContext) = DimTag[(3, gmsh.model.occ.addBox(g.corner..., g.size...))]

build!(g::Cylinder, ::BuildContext) =
    DimTag[(3, gmsh.model.occ.addCylinder(g.base..., g.axis..., g.radius))]

function build!(g::Fibre, ctx::BuildContext)
    d = g.direction
    p0 = g.point .+ d .* _dot(ctx.center .- g.point, d)   # projection of the domain centre
    H = 4 * ctx.diag
    return DimTag[(3, gmsh.model.occ.addCylinder((p0 .- H .* d)..., (2H .* d)..., g.radius))]
end

build!(g::Disk, ::BuildContext) = DimTag[(2, gmsh.model.occ.addDisk(g.center..., 0.0, g.radius, g.radius))]

function build!(g::Ellipse, ::BuildContext)
    a, b = g.radii
    # OCC requires the first radius to be the larger one
    a >= b && return DimTag[(2, gmsh.model.occ.addDisk(g.center..., 0.0, a, b))]
    t = DimTag[(2, gmsh.model.occ.addDisk(g.center..., 0.0, b, a))]
    gmsh.model.occ.rotate(t, g.center..., 0.0, 0.0, 0.0, 1.0, π / 2)
    return t
end

build!(g::Rectangle, ::BuildContext) = DimTag[(2, gmsh.model.occ.addRectangle(g.corner..., 0.0, g.size...))]

function build!(g::Fuse{D}, ctx::BuildContext) where {D}
    tags = build!(g.parts[1], ctx)
    length(g.parts) == 1 && return tags
    rest = reduce(vcat, (build!(p, ctx) for p in g.parts[2:end]))
    out, _ = gmsh.model.occ.fuse(tags, rest)
    return _nonempty(_ents(out, D), "Fuse")
end

function build!(g::Cut{D}, ctx::BuildContext) where {D}
    obj = build!(g.object, ctx)
    tools = reduce(vcat, (build!(t, ctx) for t in g.tools))
    out, _ = gmsh.model.occ.cut(obj, tools)
    return _nonempty(_ents(out, D), "Cut")
end

function build!(g::Intersect{D}, ctx::BuildContext) where {D}
    acc = build!(g.parts[1], ctx)
    for p in g.parts[2:end]
        out, _ = gmsh.model.occ.intersect(acc, build!(p, ctx))
        acc = _nonempty(_ents(out, D), "Intersect")
    end
    return acc
end

function build!(g::Translate, ctx::BuildContext)
    t = build!(g.child, ctx)
    gmsh.model.occ.translate(t, _lift3(g.shift)...)
    return t
end

function build!(g::Rotate, ctx::BuildContext)
    t = build!(g.child, ctx)
    gmsh.model.occ.rotate(t, _lift3(g.center)..., g.axis..., g.angle)
    return t
end

function build!(g::Scale{D}, ctx::BuildContext) where {D}
    t = build!(g.child, ctx)
    f = D == 3 ? g.factors : (g.factors[1], g.factors[2], 1.0)
    gmsh.model.occ.dilate(t, _lift3(g.center)..., f...)
    return t
end
