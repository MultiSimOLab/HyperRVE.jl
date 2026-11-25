const _FACE_NAMES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
_isreserved(n) = n in _FACE_NAMES || n in ("corners", "boundary") ||
                 startswith(n, "edge_") || startswith(n, "corner_")

"""
    RVE(cell, phases; matrix = :matrix, tiling = 1, overlaps = :paint)

Description of a representative volume element: a [`Cell`](@ref) filled with a background
phase called `matrix` (or nothing, `matrix = nothing`, for porous media and lattices) and the
given [`Phase`](@ref)s. Where phases overlap, the phase listed later wins
(`overlaps = :paint`); `overlaps = :error` makes any overlap an error.
"""
struct RVE{D}
    cell::Cell{D}
    phases::Vector{Phase{D}}
    matrix::Union{Nothing,String}
    tiling::NTuple{D,Int}
    overlaps::Symbol
end

function RVE(cell::Cell{D}, phases; matrix=:matrix, tiling=1, overlaps=:paint) where {D}
    D == 3 || throw(ArgumentError("2D cells are not implemented yet"))
    ps = phases isa Phase ? [phases] : collect(phases)
    all(p -> p isa Phase, ps) || throw(ArgumentError("RVE: phases must be Phase objects"))
    for p in ps
        p isa Phase{D} || throw(ArgumentError("phase $(p.name) is not $(D)D like the cell"))
    end
    m = matrix === nothing ? nothing : String(matrix)
    names = [p.name for p in ps]
    m === nothing || pushfirst!(names, m)
    for n in names
        (isempty(n) || any(isspace, n)) && throw(ArgumentError("invalid phase name $(repr(n))"))
        _isreserved(n) && throw(ArgumentError("phase name $(repr(n)) is reserved for boundary groups"))
    end
    allunique(names) || throw(ArgumentError("phase names must be unique, got $names"))
    t = tiling isa Integer ? ntuple(_ -> Int(tiling), D) :
        (length(tiling) == D && all(x -> x isa Integer, tiling) ? ntuple(i -> Int(tiling[i]), D) :
         throw(ArgumentError("tiling must be an integer or $D integers")))
    all(>=(1), t) || throw(ArgumentError("tiling must be >= 1, got $t"))
    all(==(1), t) || throw(ArgumentError("tiling (supercells) is not implemented yet"))
    overlaps in (:paint, :error) || throw(ArgumentError("overlaps must be :paint or :error, got $(repr(overlaps))"))
    for p in ps, g in p.items
        if isunbounded(g) && invariant_axis(g) == 0
            throw(ArgumentError("phase $(p.name): an unbounded geometry must be invariant along a cell " *
                                "axis (e.g. a Fibre parallel to x, y or z, possibly cut by other fibres " *
                                "parallel to it) to be a phase item. Pass crossing fibres as separate " *
                                "items, or clip the geometry with Intersect(..., Box(...)) using a box larger than the cell"))
        end
        _warn_long_cylinder(cell, p, g)
    end
    return RVE{D}(cell, Vector{Phase{D}}(ps), m, t, overlaps)
end

function _warn_long_cylinder(cell, p, g)
    g isa Cylinder || return nothing
    for d in 1:3
        if cell.periodic[d] && abs(abs(g.axis[d]) - _norm(g.axis)) < 1e-12 * _norm(g.axis) &&
           _norm(g.axis) > cell.size[d]
            @warn "phase $(p.name): a Cylinder longer than the periodic cell along axis $d overlaps " *
                  "its own periodic image; use Fibre(point, direction, radius) instead"
        end
    end
    return nothing
end

"""
    supercell(rve, n)

The same RVE repeated `n` times (integer or one integer per direction).
"""
supercell(rve::RVE{D}, n) where {D} =
    RVE(rve.cell, rve.phases; matrix=rve.matrix === nothing ? nothing : rve.matrix,
        tiling=rve.tiling .* (n isa Integer ? ntuple(_ -> Int(n), D) : Tuple(n)), overlaps=rve.overlaps)

domain_origin(rve::RVE) = rve.cell.origin
domain_size(rve::RVE) = rve.tiling .* rve.cell.size
volume(rve::RVE) = prod(domain_size(rve))

function Base.show(io::IO, r::RVE{D}) where {D}
    print(io, "RVE{$D}(", r.cell, ", phases: ",
          join((p.name for p in r.phases), ", "),
          r.matrix === nothing ? ", no matrix" : ", matrix: $(r.matrix)",
          all(==(1), r.tiling) ? "" : ", tiling $(join(r.tiling, "×"))", ")")
end
