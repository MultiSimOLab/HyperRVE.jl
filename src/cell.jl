"""
    Cell(size; origin = 0, periodic = true)
    Cell(; size, origin = 0, periodic = true)

Rectangular unit cell with edge lengths `size` (2 or 3 components) and minimum corner
`origin`. `periodic` is a `Bool` for all directions or one `Bool` per direction.
Lengths are unit-agnostic: they are written to the mesh exactly as given.

```julia
Cell((1.0, 1.0, 1.0))                                  # periodic unit cube at the origin
Cell((2.0, 1.0, 0.5); origin = (10, 0, 0), periodic = (true, true, false))
Cell((1.0, sqrt(3.0)))                                 # 2D cell
```
"""
struct Cell{D}
    origin::Vec{D}
    size::Vec{D}
    periodic::NTuple{D,Bool}
    function Cell{D}(origin::Vec{D}, size::Vec{D}, periodic::NTuple{D,Bool}) where {D}
        D in (2, 3) || throw(ArgumentError("a cell must have 2 or 3 dimensions, got $D"))
        all(>(0), size) || throw(ArgumentError("cell size must be positive, got $size"))
        return new{D}(origin, size, periodic)
    end
end

function Cell(size; origin=nothing, periodic=true)
    s = _vec_any(size, "cell size")
    D = length(s)
    o = origin === nothing ? ntuple(_ -> 0.0, D) : _vec(Val(D), origin, "cell origin")
    p = if periodic isa Bool
        ntuple(_ -> periodic, D)
    elseif periodic isa Union{Tuple,AbstractVector} && length(periodic) == D && all(x -> x isa Bool, periodic)
        ntuple(i -> periodic[i], D)
    else
        throw(ArgumentError("periodic must be a Bool or $D Bools, got $(repr(periodic))"))
    end
    return Cell{D}(o, s, p)
end
Cell(; size, origin=nothing, periodic=true) = Cell(size; origin=origin, periodic=periodic)

"""
    cellpoint(cell, ξ)

Physical coordinates of the point with fractional coordinates `ξ` (0 = minimum corner,
1 = maximum corner along each direction).
"""
cellpoint(cell::Cell{D}, ξ) where {D} = cell.origin .+ _vec(Val(D), ξ, "fractional point") .* cell.size

"""
    volume(cell)

Volume (area in 2D) of the cell.
"""
volume(c::Cell) = prod(c.size)
Base.ndims(::Cell{D}) where {D} = D

function Base.show(io::IO, c::Cell{D}) where {D}
    per = join(("xyz"[d] for d in 1:D if c.periodic[d]))
    print(io, "Cell{$D}(", join(c.size, " × "), " at ", c.origin,
          isempty(per) ? ", not periodic)" : ", periodic $per)")
end
