"""
    InterfaceRefinement(; h = nothing, level = 1.0, width = nothing, growth = nothing)
    InterfaceRefinement(h; kwargs...)

Mesh refinement on the interfaces of a phase (the surfaces it shares with other phases, not
the cell faces). The element size is `h` at the interface (if not given, `2π ℓ / (15 level)`
with `ℓ` the characteristic length of the item, the rule of the 2024 API), is kept up to the
distance `width` (default `h`) and then grows with ratio `growth` (default
`MeshOptions.growth`) until it reaches the bulk size. All lengths in the units of the cell.

In [`Phase`](@ref), `refine = 0.02` is a shorthand for `InterfaceRefinement(h = 0.02)`.
"""
struct InterfaceRefinement
    h::Union{Nothing,Float64}
    level::Float64
    width::Union{Nothing,Float64}
    growth::Union{Nothing,Float64}
    function InterfaceRefinement(h, level, width, growth)
        h = h === nothing ? nothing : _positive(h, "InterfaceRefinement h")
        if width !== nothing
            (width isa Real && isfinite(width) && width >= 0) ||
                throw(ArgumentError("InterfaceRefinement width must be >= 0, got $width"))
            width = Float64(width)
        end
        if growth !== nothing
            growth = _positive(growth, "InterfaceRefinement growth")
            growth > 1 || throw(ArgumentError("InterfaceRefinement growth must be > 1, got $growth"))
        end
        return new(h, _positive(level, "InterfaceRefinement level"), width, growth)
    end
end
InterfaceRefinement(; h=nothing, level=1.0, width=nothing, growth=nothing) =
    InterfaceRefinement(h, level, width, growth)
InterfaceRefinement(h::Real; level=1.0, width=nothing, growth=nothing) =
    InterfaceRefinement(h, level, width, growth)

_refinement(x::Nothing) = nothing
_refinement(x::Bool) = x ? InterfaceRefinement() : nothing
_refinement(x::Real) = InterfaceRefinement(h=x)
function _refinement(x::Integer)
    @warn "refine = $x is read as an element size of $x at the interfaces (in the 2024 API an " *
          "integer was a refinement level: use refine = InterfaceRefinement(level = $x) for that)"
    return InterfaceRefinement(h=x)
end
_refinement(x::InterfaceRefinement) = x
_refinement(x::NamedTuple) = InterfaceRefinement(; x...)
_refinement(x) = throw(ArgumentError("refine must be a size, true/false, a NamedTuple or an " *
                                     "InterfaceRefinement, got $(repr(x))"))

"""
    Phase(name, geometry; h = nothing, refine = nothing, fuse = false)

A named material phase. `geometry` is one [`Geometry`](@ref) or a collection of them (the
items of the phase). Items are not fused: where items of the same phase overlap (or overlap
their own periodic images) the mesh keeps the internal surfaces, which all belong to the
phase. Use [`Fuse`](@ref) to obtain one body. Each item crossing a periodic face of the cell
is completed with its periodic images.

- `h`: element size inside the phase (default: the global size of [`generate`](@ref)).
- `fuse`: unite the items of the phase and their periodic images into one body before
  meshing, so that no internal surfaces remain where items overlap (strut lattices, overlapping
  particles). Refinement by `level` then uses the smallest item of the phase.
- `refine`: refinement on the interfaces of the phase: a size (`refine = 0.02`), `true`
  (the rule of the 2024 API), a NamedTuple of [`InterfaceRefinement`](@ref) keywords or an
  `InterfaceRefinement`.
"""
struct Phase{D}
    name::String
    items::Vector{Geometry{D}}
    h::Union{Nothing,Float64}
    refine::Union{Nothing,InterfaceRefinement}
    fuse::Bool
end
function Phase(name::Union{Symbol,AbstractString}, geometry; h=nothing, refine=nothing, fuse::Bool=false)
    n = String(name)
    (isempty(n) || any(isspace, n)) && throw(ArgumentError("phase name must be non-empty and without spaces, got $(repr(n))"))
    items = geometry isa Geometry ? Geometry[geometry] : collect(geometry)
    v = _geometry_vector(items, "Phase $n")
    hp = h === nothing ? nothing : _positive(h, "Phase $n h")
    return Phase{_dim(v[1])}(n, v, hp, _refinement(refine), fuse)
end

Base.show(io::IO, p::Phase{D}) where {D} =
    print(io, "Phase{$D}(:", p.name, ", ", length(p.items), length(p.items) == 1 ? " item" : " items",
          p.h === nothing ? "" : ", h = $(p.h)", p.refine === nothing ? "" : ", refined",
          p.fuse ? ", fused" : "", ")")
