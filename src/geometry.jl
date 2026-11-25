"""
    Geometry{D}

Abstract type of all geometric descriptions in `D` dimensions. Geometries are immutable,
lazy trees: nothing is sent to Gmsh until [`generate`](@ref).

Primitives: [`Sphere`](@ref), [`Ellipsoid`](@ref), [`Box`](@ref), [`Cylinder`](@ref),
[`Fibre`](@ref). Booleans: [`Fuse`](@ref), [`Cut`](@ref), [`Intersect`](@ref) (also
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

Base.union(a::Geometry{D}, b::Geometry{D}) where {D} = Fuse(a, b)
Base.intersect(a::Geometry{D}, b::Geometry{D}) where {D} = Intersect(a, b)
Base.setdiff(a::Geometry{D}, b::Geometry{D}) where {D} = Cut(a, b)

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

Rotation by `angle` (radians) about the axis through `center`. By default a geometry rotates
about its own reference point (the centre of a sphere, ellipsoid or box, the mid-point of a
cylinder axis). In 2D the axis is always `z`.
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
    ax = D == 2 ? (0.0, 0.0, 1.0) : _normalize(_vec(Val(3), axis, "Rotate axis"), "Rotate axis")
    return Rotate{D}(g, ax, Float64(angle), _vec(Val(D), center, "Rotate center"))
end

"""
    Scale(geometry, factors; center = refpoint(geometry))

Scaling by a number or by one factor per direction about `center`.
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

# ----------------------------------------------------------------------------- traits

"""
    refpoint(g)

Reference point of a geometry: default centre of [`Rotate`](@ref) and [`Scale`](@ref).
"""
refpoint(g::Sphere) = g.center
refpoint(g::Ellipsoid) = g.center
refpoint(g::Box) = g.corner .+ g.size ./ 2
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
characteristic_length(g::Union{Cylinder,Fibre}) = g.radius
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
