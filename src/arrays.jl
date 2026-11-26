# Regular particle arrays and laminates: the classical benchmark microstructures.

const _ARRAY_SITES = Dict(
    :sc => [(0.5, 0.5, 0.5)],
    :bcc => [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
    :fcc => [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
    :square => [(0.5, 0.5)],
    :hex => [(0.0, 0.0), (0.5, 0.5)])

"""
    particle_array(cell, kind; radius = nothing, fraction = nothing) -> Vector{Sphere} or Vector{Disk}

Particles on a regular lattice in the cell: `:sc` (one sphere at the centre), `:bcc`, `:fcc`
in 3D; `:square` (one disk at the centre) and `:hex` (hexagonal packing; the cell must have the
proportions `1 : √3`) in 2D. Give the `radius` or the volume `fraction`; the particles must not
touch (fractions below about 0.52, 0.68, 0.74, 0.785 and 0.907 respectively).
`RVE(cell, Phase(:particles, particle_array(cell, :fcc; fraction = 0.3)))`.
"""
function particle_array(cell::Cell{D}, kind::Symbol; radius=nothing, fraction=nothing) where {D}
    haskey(_ARRAY_SITES, kind) && length(first(_ARRAY_SITES[kind])) == D ||
        throw(ArgumentError("particle_array: kinds are :sc, :bcc, :fcc (3D) and :square, :hex (2D), got $(repr(kind)) in $(D)D"))
    (radius === nothing) != (fraction === nothing) || throw(ArgumentError("particle_array: give radius or fraction"))
    if kind === :hex
        abs(cell.size[2] / cell.size[1] - sqrt(3)) <= 1e-6 * sqrt(3) ||
            throw(ArgumentError("particle_array: a hexagonal array needs a cell with proportions 1 : √3, e.g. Cell((a, √3 a))"))
    end
    sites = _ARRAY_SITES[kind]
    V = prod(cell.size)
    r = if radius !== nothing
        _positive(radius, "particle_array radius")
    else
        (fraction isa Real && 0 < fraction < 1) || throw(ArgumentError("particle_array: fraction must be in (0, 1)"))
        (fraction * V / (length(sites) * _ballvol(Val(D), 1.0)))^(1 / D)
    end
    # nearest-neighbour distance (physical, minimum image over the periodic lattice of sites)
    dmin = Inf
    for (i, a) in enumerate(sites), (j, b) in enumerate(sites), o in Iterators.product(ntuple(_ -> -1:1, D)...)
        i == j && all(==(0), o) && continue
        dmin = min(dmin, _norm((b .+ o .- a) .* cell.size))
    end
    2r < dmin * (1 - 1e-9) ||
        throw(ArgumentError("particle_array: particles of radius $(round(r; sigdigits=4)) touch or overlap (spacing $(round(dmin; sigdigits=4))); reduce radius or fraction"))
    return D == 3 ? Sphere[Sphere(cellpoint(cell, s), r) for s in sites] : Disk[Disk(cellpoint(cell, s), r) for s in sites]
end

"""
    laminate(cell; axis = 1, fractions, names = nothing, offset = 0.0) -> Vector{Phase}

Layers normal to the cell axis `axis` filling the cell: one phase per entry of `fractions`
(which must add up to 1), named `layer_1`, `layer_2`, … or `names`. `offset` (fraction of the
period) shifts the stack along the axis; layers crossing the face are completed periodically.
Use without matrix: `RVE(cell, laminate(cell; fractions = (0.3, 0.7)); matrix = nothing)`.
"""
function laminate(cell::Cell{D}; axis=1, fractions, names=nothing, offset=0.0) where {D}
    (axis isa Integer && !(axis isa Bool) && 1 <= axis <= D) || throw(ArgumentError("laminate: axis must be 1…$D"))
    f = collect(Float64, fractions)
    (length(f) >= 2 && all(>(0), f) && abs(sum(f) - 1) <= 1e-9) ||
        throw(ArgumentError("laminate: give at least two positive fractions adding up to 1"))
    nm = names === nothing ? [Symbol("layer_", k) for k in eachindex(f)] : collect(names)
    length(nm) == length(f) || throw(ArgumentError("laminate: one name per layer"))
    offset isa Real || throw(ArgumentError("laminate: offset must be a number"))
    L = cell.size[axis]
    z = mod(offset, 1.0)
    phases = Phase{D}[]
    for (k, fk) in enumerate(f)
        corner = ntuple(d -> d == axis ? cell.origin[d] + z * L : cell.origin[d], D)
        size = ntuple(d -> d == axis ? fk * L : cell.size[d], D)
        push!(phases, Phase(nm[k], D == 3 ? Box(corner, size) : Rectangle(corner, size)))
        z += fk
    end
    return phases
end
