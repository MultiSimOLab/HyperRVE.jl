"""
    InterfaceRefinement(; h = nothing, level = 1.0, distmin = nothing, distmax = nothing)

Mesh refinement around the interface of each item of a phase. The element size is `h` at the
interface (default `2π ℓ / (15 level)`, with `ℓ` the characteristic length of the item), is
kept up to the distance `distmin` (default `ℓ/5`) and grows to the global size at `distmax`
(default `ℓ/2.5`). All lengths in the units of the cell.
"""
struct InterfaceRefinement
    h::Union{Nothing,Float64}
    level::Float64
    distmin::Union{Nothing,Float64}
    distmax::Union{Nothing,Float64}
    function InterfaceRefinement(h, level, distmin, distmax)
        _opt(x, n) = x === nothing ? nothing : _positive(x, n)
        h = _opt(h, "InterfaceRefinement h")
        dmin = _opt(distmin, "InterfaceRefinement distmin")
        dmax = _opt(distmax, "InterfaceRefinement distmax")
        dmin !== nothing && dmax !== nothing && dmin >= dmax &&
            throw(ArgumentError("InterfaceRefinement: distmin must be smaller than distmax"))
        return new(h, _positive(level, "InterfaceRefinement level"), dmin, dmax)
    end
end
InterfaceRefinement(; h=nothing, level=1.0, distmin=nothing, distmax=nothing) =
    InterfaceRefinement(h, level, distmin, distmax)

"""
    Phase(name, geometry; refine = nothing)

A named material phase. `geometry` is one [`Geometry`](@ref) or a collection of them (the
items of the phase). Items are not fused: where items of the same phase overlap (or overlap
their own periodic images) the mesh keeps the internal surfaces, which all belong to the
phase. Use [`Fuse`](@ref) to obtain one body. Each item crossing a periodic face of the cell
is completed with its periodic images. `refine` is an [`InterfaceRefinement`](@ref) or `nothing`
(uniform global size).
"""
struct Phase{D}
    name::String
    items::Vector{Geometry{D}}
    refine::Union{Nothing,InterfaceRefinement}
end
function Phase(name::Union{Symbol,AbstractString}, geometry; refine=nothing)
    n = String(name)
    (isempty(n) || any(isspace, n)) && throw(ArgumentError("phase name must be non-empty and without spaces, got $(repr(n))"))
    items = geometry isa Geometry ? Geometry[geometry] : collect(geometry)
    v = _geometry_vector(items, "Phase $n")
    refine isa Union{Nothing,InterfaceRefinement} ||
        throw(ArgumentError("Phase $n: refine must be an InterfaceRefinement or nothing"))
    return Phase{_dim(v[1])}(n, v, refine)
end

Base.show(io::IO, p::Phase{D}) where {D} =
    print(io, "Phase{$D}(:", p.name, ", ", length(p.items), length(p.items) == 1 ? " item" : " items",
          p.refine === nothing ? "" : ", refined", ")")
