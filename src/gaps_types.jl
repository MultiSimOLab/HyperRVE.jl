# Local refinement of thin gaps between inclusions (and between an inclusion and a cell face):
# a mesh size callback that follows the local width of the ligament, evaluated in Julia.

"""
    GapRefinement(; elements = 2, faces = true, growth = nothing)

Refinement of the thin gaps between spheres or disks (any phases, void ones included) and,
with `faces = true`, between them and the cell faces: `elements` elements across each gap at
its narrowest point, the size following the local width of the ligament and growing with
`growth` (default `MeshOptions.growth`) into the inclusions. Only gaps narrower than the
smaller radius and than `elements` times the bulk size act, and only with a matrix (without
one the gaps are empty). A gap cut by a periodic face counts as two halves (`elements / 2`
between the inclusion and the face, at least 1). Use it as `MeshOptions(gaps = 2)` or
`MeshOptions(gaps = GapRefinement(elements = 3))`, with the `:delaunay` (default) or `:hxt`
3D algorithm; `optimize = :netgen` removes the occasional sliver in a ligament.
"""
struct GapRefinement
    elements::Float64
    faces::Bool
    growth::Union{Nothing,Float64}
    function GapRefinement(; elements=2, faces::Bool=true, growth=nothing)
        e = _positive(elements, "GapRefinement elements")
        g = growth === nothing ? nothing : _positive(growth, "GapRefinement growth")
        g === nothing || g > 1 || throw(ArgumentError("GapRefinement: growth must be > 1, got $g"))
        return new(Float64(e), faces, g)
    end
end
Base.show(io::IO, g::GapRefinement) =
    print(io, "GapRefinement(elements = ", g.elements, ", faces = ", g.faces,
          g.growth === nothing ? "" : ", growth = $(g.growth)", ")")
Base.:(==)(a::GapRefinement, b::GapRefinement) = a.elements == b.elements && a.faces == b.faces && a.growth == b.growth

_gaps(::Nothing) = nothing
_gaps(x::Bool) = x ? GapRefinement() : nothing
_gaps(x::GapRefinement) = x
_gaps(x::Real) = GapRefinement(elements=x)
_gaps(x::NamedTuple) = GapRefinement(; x...)
_gaps(x) = throw(ArgumentError("MeshOptions: gaps must be nothing, true/false, a number of elements, a named tuple " *
                               "or a GapRefinement, got $(repr(x))"))

