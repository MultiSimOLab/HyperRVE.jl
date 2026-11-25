"""
    MeshOptions(; h = nothing, curvature = 12, order = 1, algorithm = :default, optimize = true,
                recombine = false)

Global meshing options. `h` is the global element size in the units of the cell (default:
smallest cell edge / 10). `curvature` is the minimum number of elements per 2π along curved
boundaries, so that small inclusions are resolved even when `h` is large (0 disables it).
"""
Base.@kwdef struct MeshOptions
    h::Union{Nothing,Float64} = nothing
    curvature::Int = 12
    order::Int = 1
    algorithm::Symbol = :default
    optimize::Bool = true
    recombine::Bool = false
end

function _check(o::MeshOptions)
    o.h === nothing || _positive(o.h, "MeshOptions h")
    o.curvature >= 0 || throw(ArgumentError("MeshOptions: curvature must be >= 0"))
    o.optimize || throw(ArgumentError("MeshOptions: optimize = false is not implemented yet"))
    o.order == 1 || throw(ArgumentError("MeshOptions: order = $(o.order) is not implemented yet"))
    o.algorithm == :default || throw(ArgumentError("MeshOptions: algorithm = $(o.algorithm) is not implemented yet"))
    o.recombine && throw(ArgumentError("MeshOptions: recombine is not implemented yet"))
    return o
end
