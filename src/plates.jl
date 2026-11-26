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
_pvol(P) = length(P.faces) < 4 ? 0.0 : abs(_polyhedron_volume(_indexed(P, 1e-12)...))
function _half(P, n, c, eps)
    Q = _Poly3([copy(f) for f in P.faces], copy(P.labels))
    _clip!(Q, n, c, (0, (0, 0, 0)), eps)
    return Q
end

# convex cells left in the fractional box [a, b]^3 (cell-local coordinates) by the plates
function _plate_voids(L, k, t, a, b)
    eps = 1e-11 * maximum(L); vmin = 1e-12 * prod(L)
    pieces = [_box_poly(a .* L, b .* L)]
    for (kind, th) in zip(k, t), m in _PLATE_PLANES[kind]
        g = m ./ L; ng = _norm(g); n = g ./ ng
        s0 = sum(min(m[i] * a, m[i] * b) for i in 1:3); s1 = sum(max(m[i] * a, m[i] * b) for i in 1:3)
        for q in ceil(Int, s0 - 1e-9):floor(Int, s1 + 1e-9)
            new = _Poly3[]
            for P in pieces, Q in (_half(P, n, q / ng - th / 2, eps), _half(P, .-n, -(q / ng + th / 2), eps))
                _pvol(Q) > vmin && push!(new, Q)
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
`1 − (1 − x_sc − 2√2 x_bcc)³`; the first-order estimates `3x`, `6√2 x`, `4√3 x` overestimate.
"""
plate_density(kind, x) = plate_density(Cell((1.0, 1.0, 1.0)), kind; thickness=x)
function plate_density(cell::Cell{3}, kind; thickness)
    k = _plate_kinds(kind); t = _plate_thicknesses(k, thickness)
    return 1 - sum(_pvol, _plate_voids(cell.size, k, t, 0.0, 1.0); init=0.0) / volume(cell)
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
    k = _plate_kinds(kind); r = _plate_ratio(k, ratio)
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
or as a second material (`Phase(:filler, cells)`). Very thin plates need `h` about `t / 2`.
"""
function plate_cells(cell::Cell{3}, kind; thickness=nothing, density=nothing, ratio=:isotropic)
    (thickness === nothing) == (density === nothing) && throw(ArgumentError("plate_cells: give exactly one of thickness or density"))
    k = _plate_kinds(kind)
    density === nothing || (density isa Real && 0 < density < 1) ||
        throw(ArgumentError("plate_cells: density must be in (0, 1), got $(repr(density))"))
    t = density === nothing ? _plate_thicknesses(k, thickness) : Tuple(plate_thickness(cell, k, density; ratio=ratio))
    L = cell.size; per = cell.periodic
    eps = 1e-11 * maximum(L)
    out = Polyhedron[]
    for P in _plate_voids(L, k, t, -1.0, 2.0)
        ξ = _polyhedron_centroid(_indexed(P, 1e-12)...) ./ L
        all(d -> !per[d] || -1e-9 <= ξ[d] < 1 - 1e-9, 1:3) || continue       # one representative per periodic class
        for d in 1:3                                                        # free faces cut the cells
            per[d] && continue
            e = ntuple(i -> i == d ? 1.0 : 0.0, 3)
            _clip!(P, e, L[d], (0, (0, 0, 0)), eps); _clip!(P, .-e, 0.0, (0, (0, 0, 0)), eps)
        end
        _pvol(P) > 1e-9 * prod(L) || continue
        v, f = _indexed(P, 1e-9 * maximum(L))
        any(p -> any(x -> x < -1 + 1e-9 || x > 2 - 1e-9, p ./ L), v) && error("plate_cells: a cell is larger than the cell of the lattice")
        push!(out, Polyhedron([cell.origin .+ p for p in v], f))
    end
    isempty(out) && throw(ArgumentError("plate_cells: the plates fill the cell (thickness $t)"))
    return out
end

