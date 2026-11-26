# Plate lattices (closed-cell): plates on the simple cubic, body-centred and face-centred cubic
# planes and their combinations. The plates are the complement of the convex closed cells, which
# are built by clipping (as the Voronoi cells) and returned as polyhedra.

# plane families (Miller indices): planes m · ξ = k, k integer, ξ fractional coordinates
const _PLATE_PLANES = Dict(:sc => [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                           :bcc => [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)],
                           :fcc => [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)])
# thickness ratio of the elastically isotropic SC + FCC ("cubic + octet") plate lattice in the
# thin-plate limit, t_sc / t_fcc = 8√3/9 (Berger, Wadley and McMeeking, Nature 543, 2017)
const _PLATE_ISOTROPIC = Dict((:sc, :fcc) => (8sqrt(3) / 9, 1.0))

function _plate_kinds(kind)
    k = kind isa Symbol ? (kind,) : Tuple(kind)
    (!isempty(k) && all(x -> x isa Symbol && haskey(_PLATE_PLANES, x), k) && allunique(k)) ||
        throw(ArgumentError("plate lattice: kind must be :sc, :bcc, :fcc or a tuple of different ones, got $(repr(kind))"))
    return k
end
function _plate_ratio(k, ratio)
    length(k) == 1 && return (1.0,)
    if ratio === :isotropic
        for (key, r) in _PLATE_ISOTROPIC
            k == key && return r
            k == reverse(key) && return reverse(r)
        end
        throw(ArgumentError("plate lattice: ratio = :isotropic is known only for (:sc, :fcc); give ratio = (r1, r2, ...) or the thicknesses"))
    end
    (ratio isa Union{Tuple,AbstractVector} && length(ratio) == length(k) && all(x -> x isa Real && x > 0, ratio)) ||
        throw(ArgumentError("plate lattice: ratio must be :isotropic or $(length(k)) positive numbers"))
    return Tuple(Float64.(ratio))
end
function _plate_thicknesses(k, thickness)
    t = thickness isa Real ? ntuple(_ -> thickness, length(k)) : Tuple(thickness)
    length(t) == length(k) || throw(ArgumentError("plate lattice: thickness must be a number or one per family ($(length(k)))"))
    return map(x -> _positive(x, "plate thickness"), t)
end

function _box_poly(lo, hi)
    c(a, b, cc) = (a == 0 ? lo[1] : hi[1], b == 0 ? lo[2] : hi[2], cc == 0 ? lo[3] : hi[3])
    F = [[c(0, 0, 0), c(0, 0, 1), c(0, 1, 1), c(0, 1, 0)], [c(1, 0, 0), c(1, 1, 0), c(1, 1, 1), c(1, 0, 1)],
         [c(0, 0, 0), c(1, 0, 0), c(1, 0, 1), c(0, 0, 1)], [c(0, 1, 0), c(0, 1, 1), c(1, 1, 1), c(1, 1, 0)],
         [c(0, 0, 0), c(0, 1, 0), c(1, 1, 0), c(1, 0, 0)], [c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)]]
    return _Poly3(F, [(0, (0, 0, 0)) for _ in F])
end
# the clipping works on the cell scaled to a largest edge of 1, so the fixed tolerances are relative
_pvol(P) = length(P.faces) < 4 ? 0.0 : abs(_polyhedron_volume(_indexed(P, 1e-12)...))
function _half(P, n, c, eps)
    Q = _Poly3([copy(f) for f in P.faces], copy(P.labels))
    _clip!(Q, n, c, (0, (0, 0, 0)), eps)
    return Q
end

# runs f, turning the failures of nearly degenerate polyhedra into a clear error
function _degenerate(f)
    try
        return f()
    catch e
        e isa ArgumentError || rethrow()
        throw(ArgumentError("plate lattice: degenerate cells for these thicknesses (a plate passes almost exactly " *
                            "through the corners of the cells of another family, e.g. t_sc / t_fcc = √3); " *
                            "change the thicknesses slightly"))
    end
end

# convex cells left in the fractional box [a, b]^3 (cell-local coordinates, largest edge 1) by the plates
function _plate_voids(L, k, t, a, b)
    maximum(L) ≈ 1 || error("_plate_voids: scale the cell to a largest edge of 1")
    eps = 1e-11; vmin = 1e-10 * prod(L)
    pieces = [_box_poly(a .* L, b .* L)]
    for (kind, th) in zip(k, t), m in _PLATE_PLANES[kind]
        g = m ./ L; ng = _norm(g); n = g ./ ng
        s0 = sum(min(m[i] * a, m[i] * b) for i in 1:3); s1 = sum(max(m[i] * a, m[i] * b) for i in 1:3)
        for q in ceil(Int, s0 - 1e-9):floor(Int, s1 + 1e-9)
            new = _Poly3[]
            for P in pieces, Q in (_half(P, n, q / ng - th / 2, eps), _half(P, .-n, -(q / ng + th / 2), eps))
                _degenerate(() -> _pvol(Q)) > vmin && push!(new, Q)
            end
            pieces = new
        end
    end
    return pieces
end

"""
    plate_density(kind, x)                      # cubic cell, x = thickness / edge (one per family for a tuple)
    plate_density(cell, kind; thickness)        # any cell

Exact relative density of a plate lattice (volume of the plates, their intersections counted
once), computed by clipping. For a cube with `x = t / a`: SC `1 − (1 − x)³`, BCC
`1 − (1 − 2√2 x)³`, FCC `1 − ⅓ max(0, 1 − 2√3 x)³ − ⅔ (1 − √3 x)³`, and SC + BCC
`1 − (1 − x_sc − 2√2 x_bcc)³`, valid for `x ≤ 1` (SC), `x ≤ 1/(2√2)` (BCC), `x ≤ 1/√3` (FCC; the
tetrahedral cells vanish at `1/(2√3)`) and `x_sc + 2√2 x_bcc ≤ 1`; the first-order estimates
`3x`, `6√2 x`, `4√3 x` overestimate.
"""
plate_density(kind, x) = plate_density(Cell((1.0, 1.0, 1.0)), kind; thickness=x)
function plate_density(cell::Cell{3}, kind; thickness)
    k = _plate_kinds(kind); t = _plate_thicknesses(k, thickness)
    sc = maximum(cell.size); L = cell.size ./ sc
    return 1 - sum(_pvol, _plate_voids(L, k, t ./ sc, 0.0, 1.0); init=0.0) / prod(L)
end

"""
    plate_thickness(kind, density; ratio = :isotropic)             # cubic cell: thickness / edge
    plate_thickness(cell, kind, density; ratio = :isotropic)       # any cell: thickness

Plate thickness (one per family for a tuple of kinds, in the proportions `ratio`) that gives the
relative `density`; bisection on the exact density. `ratio = :isotropic` is the elastically
isotropic SC + FCC lattice, `t_sc / t_fcc = 8√3/9` (Berger, Wadley and McMeeking, 2017); for
other combinations give `ratio` (one number per family).
"""
plate_thickness(kind, density::Real; ratio=:isotropic) = plate_thickness(Cell((1.0, 1.0, 1.0)), kind, density; ratio=ratio)
function plate_thickness(cell::Cell{3}, kind, density::Real; ratio=:isotropic)
    0 < density < 1 || throw(ArgumentError("plate_thickness: density must be in (0, 1)"))
    k = _plate_kinds(kind)
    length(k) == 1 && ratio !== :isotropic && throw(ArgumentError("plate_thickness: ratio applies only to a tuple of kinds"))
    r = _plate_ratio(k, ratio)
    lo, hi = 0.0, minimum(cell.size) / minimum(r)
    for _ in 1:200
        m = (lo + hi) / 2
        (m == 0 ? 0.0 : plate_density(cell, k; thickness=m .* r)) < density ? (lo = m) : (hi = m)
        hi - lo < 1e-14 * maximum(cell.size) && break
    end
    t = (lo + hi) / 2 .* r
    return kind isa Symbol ? t[1] : t
end

"""
    plate_cells(cell, kind; thickness = nothing, density = nothing, ratio = :isotropic) -> Vector{Polyhedron}

The closed cells of a plate lattice: convex polyhedra, one per periodic class, left empty by
plates of thickness `thickness` (a number, or one per family) on the planes of `kind`: `:sc`
(cube faces), `:bcc` ({110}), `:fcc` ({111}) or a tuple of different ones, e.g. `(:sc, :fcc)`.
Give `thickness` or the relative `density` (with `ratio` for tuples, see
[`plate_thickness`](@ref)). Any box cell; along free directions the cells are cut by the faces.
The plates are what is left: use the cells as a void phase and the plates as the matrix,
`RVE(cell, Phase(:cells, plate_cells(cell, :fcc; density = 0.2); void = true); matrix = :plates)`,
or as a second material (`Phase(:filler, cells)`). Very thin plates need `h` about `t / 2`. Along free directions the plates on the faces
are half plates. BCC plates give lower mesh quality (min SICN around 0.01) than SC and FCC.
A warning reports cells that vanish at large thicknesses and very small features (thickness
ratios just above critical values such as `t_sc / t_fcc = √3`).
"""
function plate_cells(cell::Cell{3}, kind; thickness=nothing, density=nothing, ratio=:isotropic)
    (thickness === nothing) == (density === nothing) && throw(ArgumentError("plate_cells: give exactly one of thickness or density"))
    k = _plate_kinds(kind)
    density === nothing || (density isa Real && 0 < density < 1) ||
        throw(ArgumentError("plate_cells: density must be in (0, 1), got $(repr(density))"))
    density === nothing || ratio === :isotropic || length(k) > 1 ||
        throw(ArgumentError("plate_cells: ratio applies only to a tuple of kinds"))
    density === nothing && ratio !== :isotropic &&
        throw(ArgumentError("plate_cells: ratio applies only with density; give one thickness per family instead"))
    t = density === nothing ? _plate_thicknesses(k, thickness) : Tuple(plate_thickness(cell, k, density; ratio=ratio))
    sc = maximum(cell.size); L = cell.size ./ sc; tn = t ./ sc
    per = cell.periodic
    cells = _degenerate(() -> _plate_cells_local(L, per, k, tn))
    isempty(cells) && throw(ArgumentError("plate_cells: the plates fill the cell (thickness $t)"))
    # cells that disappear when the plates thicken (e.g. the FCC tetrahedra beyond t = a / (2√3))
    n0 = length(_degenerate(() -> _plate_cells_local(L, per, k, tn .* 1e-3)))
    length(cells) < n0 && @warn "plate_cells: $(n0 - length(cells)) of the $n0 cells of the thin-plate lattice have " *
                                "vanished at these thicknesses"
    # very small features (thickness ratios just above a critical value) give poor meshes
    emin = minimum(_min_edge(v, f) for (v, f) in cells)
    emin < 0.05 * minimum(tn) && @warn "plate_cells: the cells have edges of $(emin * sc), much shorter than the plates " *
                                      "are thick; the mesh will be poor (move the thickness ratio away from critical values such as t_sc / t_fcc = √3)"
    return _degenerate(() -> [Polyhedron([cell.origin .+ sc .* p for p in v], f) for (v, f) in cells])
end
plate_cells(::Cell{2}, kind; kwargs...) = throw(ArgumentError("plate_cells: plate lattices need a 3D cell"))

_min_edge(v, f) = minimum(_norm(v[face[i]] .- v[face[mod1(i + 1, length(face))]]) for face in f for i in eachindex(face))

# the cells (vertices, faces) in the cell scaled to a largest edge of 1
function _plate_cells_local(L, per, k, t)
    out = Tuple{Vector{NTuple{3,Float64}},Vector{Vector{Int}}}[]
    for P in _plate_voids(L, k, t, -1.0, 2.0)
        ξ = _polyhedron_centroid(_indexed(P, 1e-12)...) ./ L
        all(d -> !per[d] || -1e-9 <= ξ[d] < 1 - 1e-9, 1:3) || continue       # one representative per periodic class
        for d in 1:3                                                        # free faces cut the cells
            per[d] && continue
            e = ntuple(i -> i == d ? 1.0 : 0.0, 3)
            _clip!(P, e, L[d], (0, (0, 0, 0)), 1e-11); _clip!(P, .-e, 0.0, (0, (0, 0, 0)), 1e-11)
        end
        _pvol(P) > 1e-10 * prod(L) || continue
        v, f = _indexed(P, 1e-9)
        push!(out, ([Tuple(Float64.(p)) for p in v], f))
    end
    return out
end

