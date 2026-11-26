# Small fixed-size vector helpers on NTuple{D,Float64} (no StaticArrays dependency).

"Point or vector with `D` components."
const Vec{D} = NTuple{D,Float64}

"Gmsh `(dim, tag)` pair."
const DimTag = Tuple{Int32,Int32}

function _vec(::Val{D}, x, name::AbstractString="vector") where {D}
    x isa Union{Tuple,AbstractVector} ||
        throw(ArgumentError("$name: expected a tuple or vector of $D numbers, got $(typeof(x))"))
    length(x) == D || throw(ArgumentError("$name: expected $D components, got $(length(x))"))
    v = ntuple(i -> Float64(x[i]), Val(D))
    all(isfinite, v) || throw(ArgumentError("$name: components must be finite, got $v"))
    return v
end

function _vec_any(x, name::AbstractString="vector")
    n = x isa Union{Tuple,AbstractVector} ? length(x) : 0
    n == 2 && return _vec(Val(2), x, name)
    n == 3 && return _vec(Val(3), x, name)
    throw(ArgumentError("$name: expected 2 or 3 components, got $(repr(x))"))
end

function _positive(x, name::AbstractString)
    x isa Real || throw(ArgumentError("$name must be a real number, got $(repr(x))"))
    v = Float64(x)
    (isfinite(v) && v > 0) || throw(ArgumentError("$name must be positive and finite, got $v"))
    return v
end

_dot(a::NTuple{N}, b::NTuple{N}) where {N} = sum(a .* b)
_norm(a::Tuple) = sqrt(_dot(a, a))
_cross(a::Vec{3}, b::Vec{3}) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

function _normalize(a::Tuple, name::AbstractString="direction")
    n = _norm(a)
    n > 0 || throw(ArgumentError("$name must be a non-zero vector"))
    return a ./ n
end

_lift3(x::Vec{3}) = x
_lift3(x::Vec{2}) = (x[1], x[2], 0.0)
_drop(::Val{3}, x::Vec{3}) = x
_drop(::Val{2}, x::Vec{3}) = (x[1], x[2])

"Rotate point `p` by `angle` (radians) about the axis through `center` with unit direction `axis` (Rodrigues)."
function _rotate_point(p::Vec{D}, axis::Vec{3}, angle::Real, center::Vec{D}) where {D}
    v = _lift3(p) .- _lift3(center)
    c, s = cos(angle), sin(angle)
    r = v .* c .+ _cross(axis, v) .* s .+ axis .* (_dot(axis, v) * (1 - c))
    return _drop(Val(D), r .+ _lift3(center))
end
