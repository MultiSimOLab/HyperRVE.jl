"""
    MeshOptions(; h = nothing, hmin = nothing, hmax = nothing, growth = 1.3, curvature = 6,
                periodic_sizes = true, order = 1, algorithm = :default, optimize = true,
                recombine = false)

Global meshing options, all lengths in the units of the cell.

- `h`: element size in the bulk of the phases (default: smallest cell edge / 10). A phase
  can override it with `Phase(...; h = ...)`.
- `hmin`, `hmax`: hard lower and upper bounds of every element size.
- `growth`: ratio between the sizes of neighbouring element layers away from a refined
  interface (the size grows linearly with the distance, slope `growth - 1`).
- `curvature`: minimum number of elements per 2π along curved boundaries: near a boundary of
  radius of curvature `R` the size is at most `2πR/curvature`, so that inclusions smaller
  than `h` are still resolved (with the default 6 this only acts where `R < h`; 0 disables it).
- `periodic_sizes`: make the size field periodic (refinement near a face also refines the
  opposite face). Each refined interface is evaluated at up to 3^D - 1 shifted points, so
  meshing with refined inclusions takes about 1.5-5 times longer than without it (a larger
  `growth` shortens the transitions and reduces the cost). Only meant to be switched off for
  comparisons.
"""
struct MeshOptions
    h::Union{Nothing,Float64}
    hmin::Union{Nothing,Float64}
    hmax::Union{Nothing,Float64}
    growth::Float64
    curvature::Int
    periodic_sizes::Bool
    order::Int
    algorithm::Symbol
    optimize::Bool
    recombine::Bool
end

function MeshOptions(; h=nothing, hmin=nothing, hmax=nothing, growth=1.3, curvature=6,
                     periodic_sizes=true, order=1, algorithm=:default, optimize=true, recombine=false)
    _opt(x, n) = x === nothing ? nothing : _positive(x, n)
    h, hmin, hmax = _opt(h, "MeshOptions h"), _opt(hmin, "MeshOptions hmin"), _opt(hmax, "MeshOptions hmax")
    hmin !== nothing && hmax !== nothing && hmin > hmax &&
        throw(ArgumentError("MeshOptions: hmin = $hmin is larger than hmax = $hmax"))
    h !== nothing && hmin !== nothing && hmin > h &&
        throw(ArgumentError("MeshOptions: hmin = $hmin is larger than h = $h"))
    h !== nothing && hmax !== nothing && h > hmax &&
        throw(ArgumentError("MeshOptions: h = $h is larger than hmax = $hmax"))
    g = _positive(growth, "MeshOptions growth")
    g > 1 || throw(ArgumentError("MeshOptions: growth must be > 1, got $g"))
    curvature isa Real && isinteger(curvature) && curvature >= 0 ||
        throw(ArgumentError("MeshOptions: curvature must be a whole number >= 0, got $curvature"))
    return MeshOptions(h, hmin, hmax, g, Int(curvature), Bool(periodic_sizes), Int(order), Symbol(algorithm),
                       Bool(optimize), Bool(recombine))
end

"Copy of `o` with some fields replaced."
_with(o::MeshOptions; kw...) =
    MeshOptions(; (f => getfield(o, f) for f in fieldnames(MeshOptions))..., kw...)

function _check(o::MeshOptions)
    o.optimize || throw(ArgumentError("MeshOptions: optimize = false is not implemented yet"))
    o.order == 1 || throw(ArgumentError("MeshOptions: order = $(o.order) is not implemented yet"))
    o.algorithm == :default || throw(ArgumentError("MeshOptions: algorithm = $(o.algorithm) is not implemented yet"))
    o.recombine && throw(ArgumentError("MeshOptions: recombine is not implemented yet"))
    return o
end
