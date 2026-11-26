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
    length(sites) == 1 || all(cell.periodic) ||
        throw(ArgumentError("particle_array: $(repr(kind)) has particles on the faces and needs a periodic cell (use a supercell of a periodic cell, or :sc / :square)"))
    V = prod(cell.size)
    r = if radius !== nothing
        _positive(radius, "particle_array radius")
    else
        (fraction isa Real && 0 < fraction < 1) || throw(ArgumentError("particle_array: fraction must be in (0, 1)"))
        (fraction * V / (length(sites) * _ballvol(Val(D), 1.0)))^(1 / D)
    end
    # nearest-neighbour distance (physical, minimum image over the periodic lattice of sites)
    dmin = Inf
    for (i, a) in enumerate(sites), (j, b) in enumerate(sites), o in Iterators.product(ntuple(d -> cell.periodic[d] ? (-1:1) : (0:0), D)...)
        i == j && all(==(0), o) && continue
        dmin = min(dmin, _norm((b .+ o .- a) .* cell.size))
    end
    all(d -> cell.periodic[d] || r < cell.size[d] / 2, 1:D) ||
        throw(ArgumentError("particle_array: the particle does not fit inside the free faces of the cell"))
    2r < dmin * (1 - 1e-9) ||
        throw(ArgumentError("particle_array: particles of radius $(round(r; sigdigits=4)) touch or overlap (spacing $(round(dmin; sigdigits=4))); reduce radius or fraction"))
    return D == 3 ? Sphere[Sphere(cellpoint(cell, s), r) for s in sites] : Disk[Disk(cellpoint(cell, s), r) for s in sites]
end

"""
    laminate(cell; axis = 1, fractions, names = nothing, offset = nothing) -> Vector{Phase}

Layers normal to the cell axis `axis` filling the cell: one entry of `fractions` (adding up to
1) per layer, named `layer_1`, `layer_2`, … or `names`; layers with the same name form one
phase (`names = (:a, :b, :a, :b)`). `offset` (fraction of the period) is the position of the
start of the first layer; by default the first layer is centred on the cell face, so that no
interface lies on a face (an interface on a periodic face gets no interface group). Layers
crossing the face are split into two boxes. Use without matrix:
`RVE(cell, laminate(cell; fractions = (0.3, 0.7)); matrix = nothing)`.
"""
function laminate(cell::Cell{D}; axis=1, fractions=nothing, names=nothing, offset=nothing) where {D}
    (axis isa Integer && !(axis isa Bool) && 1 <= axis <= D) || throw(ArgumentError("laminate: axis must be 1…$D"))
    (fractions !== nothing && all(x -> x isa Real, fractions)) ||
        throw(ArgumentError("laminate: fractions must be numbers adding up to 1"))
    f = collect(Float64, fractions)
    (length(f) >= 2 && all(>(0), f) && abs(sum(f) - 1) <= 1e-9) ||
        throw(ArgumentError("laminate: give at least two positive fractions adding up to 1"))
    nm = names === nothing ? [Symbol("layer_", k) for k in eachindex(f)] :
         (names isa Union{AbstractVector,Tuple} && all(x -> x isa Union{Symbol,AbstractString}, names)) ? collect(names) :
         throw(ArgumentError("laminate: names must be a vector or tuple of symbols or strings"))
    length(nm) == length(f) || throw(ArgumentError("laminate: one name per layer"))
    offset === nothing || offset isa Real || throw(ArgumentError("laminate: offset must be a number"))
    L = cell.size[axis]
    z = offset === nothing ? -f[1] / 2 : Float64(offset)
    items = Dict{String,Vector{Geometry{D}}}()
    order = String[]
    box(z0, z1) = (corner = ntuple(d -> d == axis ? cell.origin[d] + z0 * L : cell.origin[d], D);
                   size = ntuple(d -> d == axis ? (z1 - z0) * L : cell.size[d], D);
                   D == 3 ? Box(corner, size) : Rectangle(corner, size))
    for (k, fk) in enumerate(f)
        a = mod(z, 1.0); b = a + fk
        key = String(nm[k])
        haskey(items, key) || (push!(order, key); items[key] = Geometry{D}[])
        if b <= 1 + 1e-12
            push!(items[key], box(a, min(b, 1.0)))
        else                                         # crosses the face: split into two boxes
            push!(items[key], box(a, 1.0))
            push!(items[key], box(0.0, b - 1))
        end
        z += fk
    end
    return [Phase(n, items[n]) for n in order]
end
