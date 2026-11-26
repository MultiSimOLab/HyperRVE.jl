"""
    MeshOptions(; h = nothing, hmin = nothing, hmax = nothing, growth = 1.3, curvature = nothing,
                periodic_sizes = true, order = 1, curved = true, algorithm = :default,
                surface_algorithm = :default, optimize = :auto, optimize_curved = :auto)

Global meshing options, all lengths in the units of the cell.

- `h`: element size in the bulk of the phases (default: smallest cell edge / 10). A phase
  can override it with `Phase(...; h = ...)`.
- `hmin`, `hmax`: hard lower and upper bounds of every element size.
- `growth`: ratio between the sizes of neighbouring element layers away from a refined
  interface (the size grows linearly with the distance, slope `growth - 1`).
- `curvature`: minimum number of elements per 2π along curved boundaries: near a boundary of
  radius of curvature `R` the size is at most `2πR/curvature`, so that inclusions smaller
  than `h` are still resolved (default 6, which only acts where `R < h`; 12 for curved
  order-2 meshes; 0 disables it).
- `periodic_sizes`: make the size field periodic (refinement near a face also refines the
  opposite face). Each refined interface is evaluated at up to 3^D - 1 shifted points, so
  meshing with refined inclusions takes about 1.5-5 times longer than without it (a larger
  `growth` shortens the transitions and reduces the cost). Only meant to be switched off for
  comparisons.
- `order`: 1 (linear) or 2 (quadratic elements; with `curved = true` the mid-nodes lie on the
  curved interfaces, with `curved = false` the elements are straight-sided). Order-2 meshes are
  periodic in Gmsh's sense, but GridapGmsh (0.7) reads periodic links only for first-order
  meshes: for Gridap/HyperFEM with periodic cells, generate an order-1 mesh and use
  higher-order finite element spaces on it.
- `algorithm`: 3D mesher, `:delaunay` (default), `:frontal` or `:hxt`; `surface_algorithm`:
  `:default`, `:frontal`, `:delaunay`, `:meshadapt` or `:automatic` (also used for 2D cells).
  `:default` is Frontal-Delaunay, with MeshAdapt on the B-spline surfaces (ellipsoids, bodies
  transformed with `Scale`) of inclusions inside the cell, where Frontal-Delaunay gives poor
  elements, and on any surface Frontal-Delaunay leaves without elements (Gmsh 4.15). With an
  explicit algorithm there is no such fallback; an incomplete mesh is always an error.
- `optimize`: optimisation of the linear mesh, `:default` (`true`), `:netgen` or `:none`
  (`false`); `:auto` (default) is `:netgen` for curved order-2 meshes (it removes the linear
  elements that would become inverted once curved) and `:default` otherwise. The Netgen
  optimiser acts with the `:delaunay` and `:frontal` algorithms (not with `:hxt`).
- `optimize_curved`: untangling of curved order-2 elements. `:auto` (default) is `:elastic`
  for curved order 2 and `:none` otherwise; `:elastic` removes inverted curved elements and
  keeps the periodic faces matching; `:optimize` (Gmsh's general optimiser) moves nodes of
  periodic faces, so it is only accepted for cells without periodic directions; `:none`
  leaves the curved elements as projected (coarse meshes may then contain inverted elements:
  `generate` warns and `res.stats.ninverted` counts them).
"""
struct MeshOptions
    h::Union{Nothing,Float64}
    hmin::Union{Nothing,Float64}
    hmax::Union{Nothing,Float64}
    growth::Float64
    curvature::Int
    periodic_sizes::Bool
    order::Int
    curved::Bool
    algorithm::Symbol
    surface_algorithm::Symbol
    optimize::Symbol
    optimize_curved::Symbol
end

const _ALG3D = Dict(:delaunay => 1, :frontal => 4, :hxt => 10)
const _ALG2D = Dict(:default => 6, :frontal => 6, :delaunay => 5, :meshadapt => 1, :automatic => 2)
const _HO = Dict(:optimize => "HighOrder", :elastic => "HighOrderElastic")

function MeshOptions(; h=nothing, hmin=nothing, hmax=nothing, growth=1.3, curvature=nothing,
                     periodic_sizes=true, order=1, curved=true, algorithm=:default,
                     surface_algorithm=:default, optimize=:auto, optimize_curved=:auto)
    hocurved = order == 2 && curved
    curvature === nothing && (curvature = hocurved ? 12 : 6)
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
    order in (1, 2) || throw(ArgumentError("MeshOptions: order must be 1 or 2 (higher orders are not " *
                                           "supported: periodic matching and the Gridap reader are only tested up to 2)"))
    alg = algorithm == :default ? :delaunay : Symbol(algorithm)
    haskey(_ALG3D, alg) || throw(ArgumentError("MeshOptions: algorithm must be one of $(sort(collect(keys(_ALG3D)))), got $algorithm"))
    salg = Symbol(surface_algorithm)
    haskey(_ALG2D, salg) ||
        throw(ArgumentError("MeshOptions: surface_algorithm must be one of $(sort(collect(keys(_ALG2D)))), got $surface_algorithm"))
    opt = optimize === true ? :default : optimize === false ? :none : Symbol(optimize)
    opt == :auto && (opt = hocurved ? :netgen : :default)
    opt in (:default, :netgen, :none) ||
        throw(ArgumentError("MeshOptions: optimize must be true/false or :auto, :default, :netgen, :none, got $optimize"))
    oc = Symbol(optimize_curved)
    oc == :auto && (oc = order == 2 && curved ? :elastic : :none)
    (oc == :none || haskey(_HO, oc)) ||
        throw(ArgumentError("MeshOptions: optimize_curved must be :auto, :none, :elastic or :optimize, got $optimize_curved"))
    oc != :none && !(order == 2 && curved) &&
        throw(ArgumentError("MeshOptions: optimize_curved = :$oc needs order = 2 and curved = true"))
    return MeshOptions(h, hmin, hmax, g, Int(curvature), Bool(periodic_sizes), Int(order), Bool(curved),
                       alg, salg, opt, oc)
end

"Copy of `o` with some fields replaced."
_with(o::MeshOptions; kw...) =
    MeshOptions(; (f => getfield(o, f) for f in fieldnames(MeshOptions))..., kw...)

