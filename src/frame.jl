# Coordinate frame in which the geometry is built and meshed.
# In this version the frame is the identity (physical coordinates); every tolerance is
# relative to `scale`, the largest edge of the meshed domain.

struct Frame{D}
    origin::Vec{D}     # physical coordinates of the frame origin
    scale::Float64     # reference length of the domain
    identity::Bool
end

Frame(rve::RVE{D}) where {D} = Frame{D}(ntuple(_ -> 0.0, D), maximum(domain_size(rve)), true)

to_frame(f::Frame, x::Tuple) = f.identity ? x : (x .- f.origin) ./ f.scale
to_frame(f::Frame, g::Geometry) = g
to_frame_length(f::Frame, l::Real) = f.identity ? Float64(l) : l / f.scale
from_frame_length(f::Frame, l::Real) = f.identity ? Float64(l) : l * f.scale
from_frame_measure(f::Frame{D}, v::Real, dim::Int=D) where {D} = f.identity ? Float64(v) : v * f.scale^dim

"Classification tolerance in frame units (relative 1e-6 of the domain size)."
const GEOM_TOL = 1e-6
"Node-matching tolerance, relative to the domain size."
const PERIODIC_TOL = 1e-8

geomtol(f::Frame) = f.identity ? GEOM_TOL * f.scale : GEOM_TOL
reflength(f::Frame) = f.identity ? f.scale : 1.0

restore_frame!(f::Frame) = nothing
