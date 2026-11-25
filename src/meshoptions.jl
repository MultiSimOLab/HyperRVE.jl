"""
    MeshOptions(; h = nothing, order = 1, algorithm = :default, optimize = true, recombine = false)

Global meshing options. `h` is the global element size in the units of the cell (default:
smallest cell edge / 10).
"""
Base.@kwdef struct MeshOptions
    h::Union{Nothing,Float64} = nothing
    order::Int = 1
    algorithm::Symbol = :default
    optimize::Bool = true
    recombine::Bool = false
end

function _check(o::MeshOptions)
    o.h === nothing || _positive(o.h, "MeshOptions h")
    o.order == 1 || throw(ArgumentError("MeshOptions: order = $(o.order) is not implemented yet"))
    o.algorithm == :default || throw(ArgumentError("MeshOptions: algorithm = $(o.algorithm) is not implemented yet"))
    o.recombine && throw(ArgumentError("MeshOptions: recombine is not implemented yet"))
    return o
end
