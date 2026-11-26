"""
    Geometry{D}

Abstract type of all geometric descriptions in `D` dimensions. Geometries are immutable,
lazy trees: nothing is sent to Gmsh until [`generate`](@ref).

Primitives: [`Sphere`](@ref), [`Ellipsoid`](@ref), [`Box`](@ref), [`Cylinder`](@ref),
[`Capsule`](@ref), [`Fibre`](@ref) in 3D; [`Disk`](@ref), [`Ellipse`](@ref), [`Rectangle`](@ref) in 2D. Booleans: [`Fuse`](@ref), [`Cut`](@ref), [`Intersect`](@ref) (also
`a ∪ b`, `a ∩ b`, `setdiff(a, b)`). Transforms: [`Translate`](@ref), [`Rotate`](@ref),
[`Scale`](@ref).
"""
abstract type Geometry{D} end

# ----------------------------------------------------------------------------- primitives

"""
    Sphere(center, radius)
    Sphere(; center, radius)
"""
struct Sphere <: Geometry{3}
    center::Vec{3}
    radius::Float64
    Sphere(center, radius) = new(_vec(Val(3), center, "Sphere center"), _positive(radius, "Sphere radius"))
end
Sphere(; center, radius) = Sphere(center, radius)

"""
    Ellipsoid(center, radii)
    Ellipsoid(; center, radii)

Axis-aligned ellipsoid with semi-axes `radii`. Use [`Rotate`](@ref) to orient it.
"""
struct Ellipsoid <: Geometry{3}
    center::Vec{3}
    radii::Vec{3}
    function Ellipsoid(center, radii)
        r = _vec(Val(3), radii, "Ellipsoid radii")
        all(>(0), r) || throw(ArgumentError("Ellipsoid radii must be positive, got $r"))
        return new(_vec(Val(3), center, "Ellipsoid center"), r)
    end
end
Ellipsoid(; center, radii) = Ellipsoid(center, radii)

"""
    Box(corner, size)
    Box(; corner, size)     Box(; center, size)

Axis-aligned box given by its minimum corner (or its centre) and its edge lengths.
"""
struct Box <: Geometry{3}
    corner::Vec{3}
    size::Vec{3}
    function Box(corner, size)
        s = _vec(Val(3), size, "Box size")
        all(>(0), s) || throw(ArgumentError("Box size must be positive, got $s"))
        return new(_vec(Val(3), corner, "Box corner"), s)
    end
end
function Box(; corner=nothing, center=nothing, size)
    (corner === nothing) == (center === nothing) &&
        throw(ArgumentError("Box: give exactly one of `corner` or `center`"))
    s = _vec(Val(3), size, "Box size")
    return corner !== nothing ? Box(corner, s) : Box(_vec(Val(3), center, "Box center") .- s ./ 2, s)
end

"""
    Cylinder(base, axis, radius)
    Cylinder(; base, axis, radius)
    Cylinder(; center, direction, length, radius)

Finite circular cylinder from `base` to `base .+ axis`. For a fibre that crosses the whole
cell use [`Fibre`](@ref) instead.
"""
struct Cylinder <: Geometry{3}
    base::Vec{3}
    axis::Vec{3}
    radius::Float64
    function Cylinder(base, axis, radius)
        a = _vec(Val(3), axis, "Cylinder axis")
        _norm(a) > 0 || throw(ArgumentError("Cylinder axis must be non-zero"))
        return new(_vec(Val(3), base, "Cylinder base"), a, _positive(radius, "Cylinder radius"))
    end
end
function Cylinder(; base=nothing, axis=nothing, center=nothing, direction=nothing, length=nothing, radius)
    if base !== nothing && axis !== nothing && center === nothing && direction === nothing && length === nothing
        return Cylinder(base, axis, radius)
    elseif base === nothing && axis === nothing && center !== nothing && direction !== nothing && length !== nothing
        d = _normalize(_vec(Val(3), direction, "Cylinder direction"), "Cylinder direction")
        a = d .* _positive(length, "Cylinder length")
        return Cylinder(_vec(Val(3), center, "Cylinder center") .- a ./ 2, a, radius)
    end
    throw(ArgumentError("Cylinder: give either (base, axis, radius) or (center, direction, length, radius)"))
end

"""
    Fibre(point, direction, radius)
    Fibre(; point, direction, radius)

Infinite circular cylinder through `point`. As a phase item it must be parallel to a cell
axis; its periodic images are created only across the other directions. Oblique fibres can be
used as tools inside booleans.
"""
struct Fibre <: Geometry{3}
    point::Vec{3}
    direction::Vec{3}
    radius::Float64
    Fibre(point, direction, radius) =
        new(_vec(Val(3), point, "Fibre point"),
            _normalize(_vec(Val(3), direction, "Fibre direction"), "Fibre direction"),
            _positive(radius, "Fibre radius"))
end
Fibre(; point, direction, radius) = Fibre(point, direction, radius)

"""
    Capsule(p1, p2, radius)
    Capsule(; p1, p2, radius)
    Capsule(; center, direction, length, radius)

Cylinder of axis `p1 → p2` closed by two hemispherical caps of the same radius (`length` is
the distance between the cap centres). A short fibre or a hand-made strut. For strut lattices
prefer cylinders plus one sphere per node: several
capsules ending at the same point give coincident spherical caps, a degenerate case for
OpenCASCADE booleans.
"""
struct Capsule <: Geometry{3}
    p1::Vec{3}
    p2::Vec{3}
    radius::Float64
    function Capsule(p1, p2, radius)
        a, b = _vec(Val(3), p1, "Capsule p1"), _vec(Val(3), p2, "Capsule p2")
        a == b && throw(ArgumentError("Capsule: p1 and p2 coincide; use Sphere"))
        return new(a, b, _positive(radius, "Capsule radius"))
    end
end
function Capsule(; p1=nothing, p2=nothing, center=nothing, direction=nothing, length=nothing, radius)
    if p1 !== nothing && p2 !== nothing && center === nothing && direction === nothing && length === nothing
        return Capsule(p1, p2, radius)
    elseif p1 === nothing && p2 === nothing && center !== nothing && direction !== nothing && length !== nothing
        d = _normalize(_vec(Val(3), direction, "Capsule direction"), "Capsule direction") .* (_positive(length, "Capsule length") / 2)
        c = _vec(Val(3), center, "Capsule center")
        return Capsule(c .- d, c .+ d, radius)
    end
    throw(ArgumentError("Capsule: give either (p1, p2, radius) or (center, direction, length, radius)"))
end

# ----------------------------------------------------------------------------- 2D primitives

"""
    Disk(center, radius)
    Disk(; center, radius)

Circular inclusion in a 2D cell.
"""
struct Disk <: Geometry{2}
    center::Vec{2}
    radius::Float64
    Disk(center, radius) = new(_vec(Val(2), center, "Disk center"), _positive(radius, "Disk radius"))
end
Disk(; center, radius) = Disk(center, radius)

"""
    Ellipse(center, radii)
    Ellipse(; center, radii)

Axis-aligned ellipse with semi-axes `radii` in a 2D cell. Use [`Rotate`](@ref) to orient it.
"""
struct Ellipse <: Geometry{2}
    center::Vec{2}
    radii::Vec{2}
    function Ellipse(center, radii)
        r = _vec(Val(2), radii, "Ellipse radii")
        all(>(0), r) || throw(ArgumentError("Ellipse radii must be positive, got $r"))
        return new(_vec(Val(2), center, "Ellipse center"), r)
    end
end
Ellipse(; center, radii) = Ellipse(center, radii)

"""
    Rectangle(corner, size)
    Rectangle(; corner, size)     Rectangle(; center, size)

Axis-aligned rectangle in a 2D cell, given by its minimum corner (or its centre) and sizes.
"""
struct Rectangle <: Geometry{2}
    corner::Vec{2}
    size::Vec{2}
    function Rectangle(corner, size)
        s = _vec(Val(2), size, "Rectangle size")
        all(>(0), s) || throw(ArgumentError("Rectangle size must be positive, got $s"))
        return new(_vec(Val(2), corner, "Rectangle corner"), s)
    end
end
function Rectangle(; corner=nothing, center=nothing, size)
    (corner === nothing) == (center === nothing) &&
        throw(ArgumentError("Rectangle: give exactly one of `corner` or `center`"))
    s = _vec(Val(2), size, "Rectangle size")
    return corner !== nothing ? Rectangle(corner, s) :
           Rectangle(_vec(Val(2), center, "Rectangle center") .- s ./ 2, s)
end

# ----------------------------------------------------------------------------- booleans

function _geometry_vector(parts, what)
    v = collect(parts)
    isempty(v) && throw(ArgumentError("$what needs at least one geometry"))
    all(p -> p isa Geometry, v) || throw(ArgumentError("$what: all parts must be geometries"))
    D = _dim(first(v))
    all(p -> _dim(p) == D, v) || throw(ArgumentError("$what: cannot mix 2D and 3D geometries"))
    return Vector{Geometry{D}}(v)
end
_dim(::Geometry{D}) where {D} = D

"""
    Fuse(parts...)     Fuse(parts::AbstractVector)      a ∪ b

Union of geometries (one body).
"""
struct Fuse{D} <: Geometry{D}
    parts::Vector{Geometry{D}}
    Fuse{D}(parts::Vector{Geometry{D}}) where {D} = new{D}(parts)
end
Fuse(parts::Geometry...) = (v = _geometry_vector(parts, "Fuse"); Fuse{_dim(v[1])}(v))
Fuse(parts::AbstractVector) = (v = _geometry_vector(parts, "Fuse"); Fuse{_dim(v[1])}(v))

"""
    Cut(object, tools...)      setdiff(object, tool)

`object` minus the union of `tools`.
"""
struct Cut{D} <: Geometry{D}
    object::Geometry{D}
    tools::Vector{Geometry{D}}
    Cut{D}(object::Geometry{D}, tools::Vector{Geometry{D}}) where {D} = new{D}(object, tools)
end
function Cut(object::Geometry{D}, tools::Geometry...) where {D}
    v = _geometry_vector(tools, "Cut tools")
    _dim(v[1]) == D || throw(ArgumentError("Cut: cannot mix 2D and 3D geometries"))
    return Cut{D}(object, v)
end

"""
    Intersect(parts...)     Intersect(parts::AbstractVector)      a ∩ b

Intersection of geometries.
"""
struct Intersect{D} <: Geometry{D}
    parts::Vector{Geometry{D}}
    Intersect{D}(parts::Vector{Geometry{D}}) where {D} = new{D}(parts)
end
Intersect(parts::Geometry...) = (v = _geometry_vector(parts, "Intersect"); Intersect{_dim(v[1])}(v))
Intersect(parts::AbstractVector) = (v = _geometry_vector(parts, "Intersect"); Intersect{_dim(v[1])}(v))

# untyped in D so that mixing 2D and 3D gives the ArgumentError of the constructors
Base.union(a::Geometry, b::Geometry, c::Geometry...) = Fuse(a, b, c...)
Base.intersect(a::Geometry, b::Geometry, c::Geometry...) = Intersect(a, b, c...)
Base.setdiff(a::Geometry, b::Geometry) = Cut(a, b)

# ----------------------------------------------------------------------------- transforms

"""
    Translate(geometry, shift)
"""
struct Translate{D} <: Geometry{D}
    child::Geometry{D}
    shift::Vec{D}
    Translate{D}(child::Geometry{D}, shift::Vec{D}) where {D} = new{D}(child, shift)
end
Translate(g::Geometry{D}, shift) where {D} = Translate{D}(g, _vec(Val(D), shift, "Translate shift"))

"""
    Rotate(geometry; angle, axis = (0, 0, 1), center = refpoint(geometry))

Rotation by `angle` (radians, right-handed about `axis`) about the axis through `center`. By
default a geometry rotates about its own reference point: the centre of a sphere, ellipsoid,
disk or box, the mid-point of a cylinder axis, the point of a fibre, the first part of a
boolean. In 2D the rotation is about `z` (counter-clockwise for positive angles); `axis` may
only be `(0, 0, ±1)` there.

```julia
Rotate(Ellipsoid((0.5, 0.5, 0.5), (0.3, 0.1, 0.1)); axis = (1, 1, 0), angle = π/4)
Rotate(Ellipse((0.5, 0.5), (0.3, 0.1)); angle = π/6)                # 2D
```
"""
struct Rotate{D} <: Geometry{D}
    child::Geometry{D}
    axis::Vec{3}
    angle::Float64
    center::Vec{D}
    Rotate{D}(child::Geometry{D}, axis::Vec{3}, angle::Float64, center::Vec{D}) where {D} =
        new{D}(child, axis, angle, center)
end
function Rotate(g::Geometry{D}; angle, axis=(0, 0, 1), center=refpoint(g)) where {D}
    angle isa Real && isfinite(angle) || throw(ArgumentError("Rotate: angle must be a finite number"))
    ax = _normalize(_vec(Val(3), axis, "Rotate axis"), "Rotate axis")
    θ = Float64(angle)
    if D == 2
        abs(abs(ax[3]) - 1) < _AXIS_TOL ||
            throw(ArgumentError("Rotate: in 2D the axis must be (0, 0, ±1), got $axis"))
        θ = ax[3] > 0 ? θ : -θ
        ax = (0.0, 0.0, 1.0)
    end
    return Rotate{D}(g, ax, θ, _vec(Val(D), center, "Rotate center"))
end

"""
    Scale(geometry, factors; center = refpoint(geometry))

Scaling by a number or by one factor per direction about `center` (by default the reference
point of the geometry, see [`Rotate`](@ref)).
"""
struct Scale{D} <: Geometry{D}
    child::Geometry{D}
    factors::Vec{D}
    center::Vec{D}
    Scale{D}(child::Geometry{D}, factors::Vec{D}, center::Vec{D}) where {D} = new{D}(child, factors, center)
end
function Scale(g::Geometry{D}, factors; center=refpoint(g)) where {D}
    f = factors isa Real ? ntuple(_ -> Float64(factors), D) : _vec(Val(D), factors, "Scale factors")
    all(x -> isfinite(x) && x != 0, f) || throw(ArgumentError("Scale factors must be finite and non-zero, got $f"))
    return Scale{D}(g, f, _vec(Val(D), center, "Scale center"))
end


"""
    Polyhedron(vertices, faces)
    Polyhedron(; vertices, faces)

Solid bounded by planar polygonal faces: `faces` are loops of indices into `vertices` (either
orientation). The surface must be closed (every edge in exactly two faces) and each face
planar (split non-planar faces into triangles). The solid need not be convex. A Voronoi grain,
a prism, a faceted particle.
"""
struct Polyhedron <: Geometry{3}
    vertices::Vector{Vec{3}}
    faces::Vector{Vector{Int}}
    function Polyhedron(vertices, faces)
        v = Vec{3}[_vec(Val(3), p, "Polyhedron vertex") for p in vertices]
        f = Vector{Int}[]
        for face in faces
            all(x -> x isa Integer, face) || throw(ArgumentError("Polyhedron: a face is a list of vertex indices, got $face"))
            push!(f, Int.(collect(face)))
        end
        length(v) >= 4 && length(f) >= 4 || throw(ArgumentError("Polyhedron: needs at least 4 vertices and 4 faces"))
        L = maximum(p -> _norm(p .- v[1]), v)
        edges = Dict{Tuple{Int,Int},Int}()
        for face in f
            length(unique(face)) == length(face) >= 3 ||
                throw(ArgumentError("Polyhedron: face $face needs at least 3 distinct vertices"))
            all(i -> 1 <= i <= length(v), face) || throw(ArgumentError("Polyhedron: face $face refers to a missing vertex"))
            n = _face_normal(v, face)
            _norm(n) > 1e-12 * L^2 || throw(ArgumentError("Polyhedron: face $face has zero area"))
            u = n ./ _norm(n)
            all(i -> abs(sum((v[i] .- v[face[1]]) .* u)) <= 1e-9 * L, face) ||
                throw(ArgumentError("Polyhedron: face $face is not planar; split it into triangles"))
            for (a, b) in zip(face, circshift(face, -1))
                k = minmax(a, b)
                edges[k] = get(edges, k, 0) + 1
            end
        end
        bad = [k for (k, c) in edges if c != 2]
        isempty(bad) || throw(ArgumentError("Polyhedron: the surface is not closed (edges $(first(bad, 3)) are not in exactly two faces)"))
        vol = abs(_polyhedron_volume(v, f))
        vol > 1e-12 * L^3 || throw(ArgumentError("Polyhedron: zero volume"))
        return new(v, f)
    end
end
Polyhedron(; vertices, faces) = Polyhedron(vertices, faces)

# Newell normal (twice the area vector) of a planar loop
function _face_normal(v, face)
    n = (0.0, 0.0, 0.0)
    for (a, b) in zip(face, circshift(face, -1))
        p, q = v[a], v[b]
        n = n .+ ((p[2] - q[2]) * (p[3] + q[3]), (p[3] - q[3]) * (p[1] + q[1]), (p[1] - q[1]) * (p[2] + q[2]))
    end
    return n
end

# signed volume and volume centroid with faces oriented consistently; orientation fixed by
# propagation over shared edges
function _oriented_faces(v, f)
    out = [copy(x) for x in f]
    nf = length(out)
    done = falses(nf); done[1] = true
    queue = [1]
    while !isempty(queue)
        i = popfirst!(queue)
        de = Set(zip(out[i], circshift(out[i], -1)))
        for j in 1:nf
            done[j] && continue
            ej = collect(zip(out[j], circshift(out[j], -1)))
            shares_same = any(e -> e in de, ej)
            shares_rev = any(e -> (e[2], e[1]) in de, ej)
            (shares_same || shares_rev) || continue
            shares_same && reverse!(out[j])
            done[j] = true; push!(queue, j)
        end
    end
    all(done) || throw(ArgumentError("Polyhedron: the faces do not form one connected surface"))
    return out
end

function _polyhedron_volume(v, f)
    of = _oriented_faces(v, f)
    vol = 0.0
    for face in of, k in 2:length(face)-1
        a, b, c = v[face[1]], v[face[k]], v[face[k+1]]
        vol += (a[1] * (b[2] * c[3] - b[3] * c[2]) - a[2] * (b[1] * c[3] - b[3] * c[1]) + a[3] * (b[1] * c[2] - b[2] * c[1])) / 6
    end
    return vol
end

function _polyhedron_centroid(v, f)
    of = _oriented_faces(v, f)
    vol = 0.0; c = (0.0, 0.0, 0.0)
    for face in of, k in 2:length(face)-1
        a, b, d = v[face[1]], v[face[k]], v[face[k+1]]
        t = (a[1] * (b[2] * d[3] - b[3] * d[2]) - a[2] * (b[1] * d[3] - b[3] * d[1]) + a[3] * (b[1] * d[2] - b[2] * d[1])) / 6
        vol += t; c = c .+ t .* (a .+ b .+ d) ./ 4
    end
    return c ./ vol
end

"""
    Polygon(vertices)
    Polygon(; vertices)

Simple polygon in a 2D cell (either orientation; self-intersection is not checked).
"""
struct Polygon <: Geometry{2}
    vertices::Vector{Vec{2}}
    function Polygon(vertices)
        v = Vec{2}[_vec(Val(2), p, "Polygon vertex") for p in vertices]
        length(v) >= 3 || throw(ArgumentError("Polygon: needs at least 3 vertices"))
        L = maximum(p -> _norm(p .- v[1]), v)
        abs(_polygon_area(v)) > 1e-12 * L^2 || throw(ArgumentError("Polygon: zero area"))
        all(i -> _norm(v[i] .- v[mod1(i + 1, length(v))]) > 1e-12 * L, eachindex(v)) ||
            throw(ArgumentError("Polygon: repeated consecutive vertices"))
        return new(v)
    end
end
Polygon(; vertices) = Polygon(vertices)

_polygon_area(v) = sum(v[i][1] * v[mod1(i + 1, length(v))][2] - v[mod1(i + 1, length(v))][1] * v[i][2] for i in eachindex(v)) / 2
function _polygon_centroid(v)
    A = _polygon_area(v); cx = 0.0; cy = 0.0
    for i in eachindex(v)
        p, q = v[i], v[mod1(i + 1, length(v))]
        w = p[1] * q[2] - q[1] * p[2]
        cx += (p[1] + q[1]) * w; cy += (p[2] + q[2]) * w
    end
    return (cx / (6A), cy / (6A))
end


# ----------------------------------------------------------------------------- traits

"""
    refpoint(g)

Reference point of a geometry: default centre of [`Rotate`](@ref) and [`Scale`](@ref).
"""
refpoint(g::Sphere) = g.center
refpoint(g::Ellipsoid) = g.center
refpoint(g::Box) = g.corner .+ g.size ./ 2
refpoint(g::Union{Disk,Ellipse}) = g.center
refpoint(g::Capsule) = (g.p1 .+ g.p2) ./ 2
refpoint(g::Rectangle) = g.corner .+ g.size ./ 2
refpoint(g::Polyhedron) = _polyhedron_centroid(g.vertices, g.faces)
refpoint(g::Polygon) = _polygon_centroid(g.vertices)
refpoint(g::Cylinder) = g.base .+ g.axis ./ 2
refpoint(g::Fibre) = g.point
refpoint(g::Union{Fuse,Intersect}) = refpoint(first(g.parts))
refpoint(g::Cut) = refpoint(g.object)
refpoint(g::Translate) = refpoint(g.child) .+ g.shift
refpoint(g::Rotate) = _rotate_point(refpoint(g.child), g.axis, g.angle, g.center)
refpoint(g::Scale) = g.center .+ g.factors .* (refpoint(g.child) .- g.center)

"""
    characteristic_length(g)

Length scale used by the default interface refinement (radius, smallest semi-axis or edge).
"""
characteristic_length(g::Sphere) = g.radius
characteristic_length(g::Ellipsoid) = minimum(g.radii)
characteristic_length(g::Box) = minimum(g.size)
characteristic_length(g::Disk) = g.radius
characteristic_length(g::Ellipse) = minimum(g.radii)
characteristic_length(g::Rectangle) = minimum(g.size)
characteristic_length(g::Polyhedron) =     # 3V / A: the inradius of a regular polyhedron
    3abs(_polyhedron_volume(g.vertices, g.faces)) / sum(f -> _norm(_face_normal(g.vertices, f)) / 2, g.faces)
characteristic_length(g::Polygon) =        # 2A / perimeter
    2abs(_polygon_area(g.vertices)) / sum(i -> _norm(g.vertices[i] .- g.vertices[mod1(i + 1, length(g.vertices))]), eachindex(g.vertices))
characteristic_length(g::Union{Cylinder,Fibre,Capsule}) = g.radius
characteristic_length(g::Union{Fuse,Intersect}) = minimum(characteristic_length, g.parts)
characteristic_length(g::Cut) = min(characteristic_length(g.object), minimum(characteristic_length, g.tools))
characteristic_length(g::Union{Translate,Rotate}) = characteristic_length(g.child)
characteristic_length(g::Scale) = characteristic_length(g.child) * minimum(abs, g.factors)

"""
    isunbounded(g)

Whether the geometry extends to infinity (contains a [`Fibre`](@ref) that is not clipped).
"""
isunbounded(::Geometry) = false
isunbounded(::Fibre) = true
isunbounded(g::Fuse) = any(isunbounded, g.parts)
isunbounded(g::Intersect) = all(isunbounded, g.parts)
isunbounded(g::Cut) = isunbounded(g.object)
isunbounded(g::Union{Translate,Rotate,Scale}) = isunbounded(g.child)

const _AXIS_TOL = 1e-12

"""
    invariant_axis(g)

`d` if `g` is invariant under translations along the cell axis `d` (an aligned fibre),
otherwise `0`.
"""
invariant_axis(::Geometry) = 0
function invariant_axis(g::Fibre)
    for d in 1:3
        abs(abs(g.direction[d]) - 1) < _AXIS_TOL && return d
    end
    return 0
end
function _common_axis(parts)
    a = invariant_axis(first(parts))
    return all(p -> invariant_axis(p) == a, parts) ? a : 0
end
invariant_axis(g::Fuse) = _common_axis(g.parts)
invariant_axis(g::Intersect) = _common_axis(g.parts)
function invariant_axis(g::Cut)
    a = invariant_axis(g.object)
    return all(t -> invariant_axis(t) == a, g.tools) ? a : 0
end
invariant_axis(g::Union{Translate,Scale}) = invariant_axis(g.child)
function invariant_axis(g::Rotate{D}) where {D}
    a = invariant_axis(g.child)
    a == 0 && return 0
    e = ntuple(i -> i == a ? 1.0 : 0.0, 3)
    return abs(abs(_dot(e, g.axis)) - 1) < _AXIS_TOL || g.angle == 0 ? a : 0
end
