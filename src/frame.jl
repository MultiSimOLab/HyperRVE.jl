# Coordinate frame in which the geometry is built and meshed.
#
# OCC works with an absolute tolerance (~1e-7) and Gmsh with Geometry.Tolerance = 1e-8, so a
# cell of 1e-6 (µm in metres) would be only a few tolerances wide. The geometry is therefore
# built and meshed in a normalised frame (domain origin at 0, longest domain edge = 1) and
# the mesh nodes are mapped back to physical coordinates before anything is exported. The
# same OCC/Gmsh computation runs whatever the physical scale.

struct Frame{D}
    origin::Vec{D}     # physical coordinates of the domain's minimum corner
    scale::Float64     # longest edge of the meshed domain
end

Frame(rve::RVE{D}) where {D} = Frame{D}(domain_origin(rve), maximum(domain_size(rve)))

_isidentity(f::Frame) = f.scale == 1 && all(iszero, f.origin)

to_frame(f::Frame, x::Tuple) = (x .- f.origin) ./ f.scale
to_frame_length(f::Frame, l::Real) = l / f.scale
from_frame_length(f::Frame, l::Real) = l * f.scale
from_frame_measure(f::Frame{D}, v::Real, dim::Int=D) where {D} = v * f.scale^dim

"Classification tolerance (frame units)."
const GEOM_TOL = 1e-6
"Node-matching tolerance, relative to the domain size (B-spline intersections are only
reproducible to ~1e-8 between a face and its translate)."
const PERIODIC_TOL = 1e-6

geomtol(::Frame) = GEOM_TOL

# Geometry in frame coordinates. A uniform scaling plus a translation maps every primitive to
# a primitive of the same kind and conjugates the transforms, so no OCC dilation is needed.
to_frame(f::Frame, g::Geometry) = _isidentity(f) ? g : _tf(f, g)

_tf(f, g::Sphere) = Sphere(to_frame(f, g.center), g.radius / f.scale)
_tf(f, g::Ellipsoid) = Ellipsoid(to_frame(f, g.center), g.radii ./ f.scale)
_tf(f, g::Box) = Box(to_frame(f, g.corner), g.size ./ f.scale)
_tf(f, g::Cylinder) = Cylinder(to_frame(f, g.base), g.axis ./ f.scale, g.radius / f.scale)
_tf(f, g::Fibre) = Fibre(to_frame(f, g.point), g.direction, g.radius / f.scale)
_tf(f, g::Fuse{D}) where {D} = Fuse{D}(Geometry{D}[_tf(f, p) for p in g.parts])
_tf(f, g::Intersect{D}) where {D} = Intersect{D}(Geometry{D}[_tf(f, p) for p in g.parts])
_tf(f, g::Cut{D}) where {D} = Cut{D}(_tf(f, g.object), Geometry{D}[_tf(f, t) for t in g.tools])
_tf(f, g::Translate{D}) where {D} = Translate{D}(_tf(f, g.child), g.shift ./ f.scale)
_tf(f, g::Rotate{D}) where {D} = Rotate{D}(_tf(f, g.child), g.axis, g.angle, to_frame(f, g.center))
_tf(f, g::Scale{D}) where {D} = Scale{D}(_tf(f, g.child), g.factors, to_frame(f, g.center))

"Map the mesh nodes from the frame back to physical coordinates."
function restore_frame!(f::Frame{D}) where {D}
    _isidentity(f) && return nothing
    s = f.scale
    o = _lift3(f.origin)
    gmsh.model.mesh.affineTransform([s, 0, 0, o[1], 0, s, 0, o[2], 0, 0, s, o[3]])
    return nothing
end
