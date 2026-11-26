"""Volume and identification of one phase in a generated RVE."""
struct PhaseInfo
    name::String
    tag::Int
    entities::Vector{Int}
    volume::Float64
    fraction::Float64
end

"""
An interface between two phases (`phases[2] == "void"` for free surfaces without matrix):
its group name and tag (0 if groups were not created), its `(D-1)` entities and its area
(length in 2D) in the units of the cell.
"""
struct InterfaceInfo
    name::String
    phases::Tuple{String,String}
    tag::Int
    entities::Vector{Int}
    area::Float64
end

"""A physical group of the generated model."""
struct GroupInfo
    name::String
    dim::Int
    tag::Int
    entities::Vector{Int}
end

"""
Summary statistics of the generated mesh (lengths in the units of the cell). Quality is the
signed inverse condition number (SICN, 1 = regular element, <= 0 = inverted);
`quality_histogram[i]` counts elements with SICN in `[(i-1)/10, i/10)` (non-positive values
in the first bin). `min_jacobian` is the minimum scaled Jacobian of order-2 elements (NaN for
order 1). `worst_point`/`worst_phase` locate the element of lowest quality and
`phase_quality` gives the minimum quality per phase. See [`quality_report`](@ref).
"""
struct MeshStats
    nnodes::Int
    nelements::Int
    order::Int
    min_quality::Float64
    mean_quality::Float64
    hmin::Float64
    hmax::Float64
    quality_histogram::NTuple{10,Int}
    ninverted::Int
    min_jacobian::Float64
    worst_point::NTuple{3,Float64}
    worst_phase::String
    phase_quality::Vector{Pair{String,Float64}}
end

"""Result of the node-by-node check of periodicity along one direction."""
struct PeriodicityCheck
    direction::Int
    npairs::Int
    max_error::Float64
    ok::Bool
end

"""
    PeriodicityError(direction, msg)

Thrown when opposite faces of a periodic cell cannot be paired or their meshes do not match.
"""
struct PeriodicityError <: Exception
    direction::Int
    msg::String
end
Base.showerror(io::IO, e::PeriodicityError) =
    print(io, "PeriodicityError along ", "xyz"[e.direction], ": ", e.msg)

"""
    RVEMesh

Result of [`generate`](@ref): phase volumes and fractions, physical groups, mesh statistics,
periodicity checks and written files. It is plain Julia data, valid after Gmsh is closed.
Query it with [`volume`](@ref), [`volume_fraction`](@ref), [`volume_fractions`](@ref),
[`phase_names`](@ref), [`group_names`](@ref), [`group_tag`](@ref) or `res[:name]`.
"""
struct RVEMesh{D}
    rve::RVE{D}
    options::MeshOptions
    phases::Vector{PhaseInfo}
    groups::Vector{GroupInfo}
    interfaces::Vector{InterfaceInfo}
    stats::Union{Nothing,MeshStats}
    periodicity::Vector{PeriodicityCheck}
    files::Vector{String}
    model::String
end

"""
    quality_report([io = stdout], res)

Print the quality statistics of a generated mesh: summary, histogram of the element quality
(SICN), minimum per phase and the location of the worst element. Works after Gmsh is closed.
"""
quality_report(r::RVEMesh) = quality_report(stdout, r)
function quality_report(io::IO, r::RVEMesh{D}) where {D}
    s = r.stats
    if s === nothing
        println(io, "no mesh (geometry_only = true)")
        return nothing
    end
    println(io, @sprintf("%d %s of order %d, %d nodes; quality (SICN) min %.3f, mean %.3f%s",
                         s.nelements, D == 3 ? "tetrahedra" : "triangles", s.order, s.nnodes,
                         s.min_quality, s.mean_quality,
                         s.order == 2 ? @sprintf(", min scaled Jacobian %.3f", s.min_jacobian) : ""))
    s.ninverted > 0 && println(io, "  ", s.ninverted, " inverted elements")
    n = max(s.nelements, 1)
    for (i, c) in enumerate(s.quality_histogram)
        bar = repeat("█", round(Int, 40 * c / n))
        println(io, @sprintf("  [%.1f, %.1f%s  %8d  %5.1f%%  %s", (i - 1) / 10, i / 10, i == 10 ? "]" : ")", c, 100c / n, bar))
    end
    for (p, q) in s.phase_quality
        println(io, @sprintf("  %-14s min quality %.3f", p, q))
    end
    println(io, @sprintf("  worst element in phase %s at (%.4g, %.4g, %.4g)", s.worst_phase, s.worst_point...))
    return nothing
end

"""    interface_names(res)  Names of the interfaces `interface_<a>_<b>` between phases."""
interface_names(r::RVEMesh) = [i.name for i in r.interfaces]

function _interface(r::RVEMesh, a, b)
    x, y = String(a), String(b)
    for i in r.interfaces
        (i.phases == (x, y) || i.phases == (y, x)) && return i
    end
    throw(KeyError("no interface between $x and $y (existing: $(join(interface_names(r), ", ")))"))
end

"""    interface_area(res, a, b)  Area (length in 2D) of the interface between phases `a` and `b` (either order; `:void` for free surfaces)."""
interface_area(r::RVEMesh, a, b) = _interface(r, a, b).area

"""    interface_tag(res, a, b)  Physical tag of the interface group between phases `a` and `b` (same in every dimension)."""
function interface_tag(r::RVEMesh, a, b)
    i = _interface(r, a, b)
    i.tag == 0 && error("interface groups were not created (generate(...; interfaces = false))")
    return i.tag
end

"""    phase_names(res)  Names of the phases (matrix first, if any)."""
phase_names(r::RVEMesh) = [p.name for p in r.phases]

function Base.getindex(r::RVEMesh, name::Union{Symbol,AbstractString})
    n = String(name)
    for p in r.phases
        p.name == n && return p
    end
    throw(KeyError("$n (available phases: $(join(phase_names(r), ", ")))"))
end

"""
    volume(res, name)   volume(res)

Volume (area in 2D) of a phase, or of the whole meshed domain, in the units of the cell,
computed by OpenCASCADE from the exact geometry. For pieces bounded by B-spline surfaces
(ellipsoids and non-uniformly scaled bodies cut by a cell face) the OCC integration is only
accurate to about 1e-2 relative; analytic bodies are accurate to round-off.
"""
volume(r::RVEMesh, name::Union{Symbol,AbstractString}) = r[name].volume
volume(r::RVEMesh) = volume(r.rve)

"""    volume_fraction(res, name)  Volume of a phase divided by the volume of the domain."""
volume_fraction(r::RVEMesh, name::Union{Symbol,AbstractString}) = r[name].fraction

"""    volume_fractions(res)  `name => fraction` for every phase."""
volume_fractions(r::RVEMesh) = [p.name => p.fraction for p in r.phases]

"""    group_names(res; dim = nothing)  Names of the physical groups (optionally of one dimension)."""
group_names(r::RVEMesh; dim=nothing) =
    unique!([g.name for g in r.groups if dim === nothing || g.dim == dim])

"""    group_tag(res, name; dim)  Physical tag of a group."""
function group_tag(r::RVEMesh, name::Union{Symbol,AbstractString}; dim::Integer)
    n = String(name)
    for g in r.groups
        g.name == n && g.dim == dim && return g.tag
    end
    throw(KeyError("$n in dimension $dim (available: $(join(group_names(r; dim=dim), ", ")))"))
end

Base.show(io::IO, r::RVEMesh{D}) where {D} =
    print(io, "RVEMesh{$D}(", length(r.phases), " phases",
          r.stats === nothing ? ", geometry only" : ", $(r.stats.nelements) elements", ")")

function Base.show(io::IO, ::MIME"text/plain", r::RVEMesh{D}) where {D}
    c = r.rve.cell
    per = join(("xyz"[d] for d in 1:D if c.periodic[d]))
    println(io, "RVEMesh{$D}: cell ", join(c.size, " × "), " at ", c.origin,
            all(==(1), r.rve.tiling) ? "" : ", tiling $(join(r.rve.tiling, "×"))",
            isempty(per) ? ", not periodic" : ", periodic $per")
    println(io, @sprintf("  %-14s %4s  %-12s  %s", "phase", "tag", "volume", "fraction"))
    for p in r.phases
        println(io, @sprintf("  %-14s %4d  %-12.5e  %.4f", p.name, p.tag, p.volume, p.fraction))
    end
    if r.stats === nothing
        println(io, "  mesh: not generated (geometry only)")
    else
        s = r.stats
        println(io, @sprintf("  mesh: %d nodes, %d %s (order %d), min quality %.3f, h ∈ [%.3g, %.3g]",
                             s.nnodes, s.nelements, D == 3 ? "tetrahedra" : "triangles", s.order,
                             s.min_quality, s.hmin, s.hmax),
                s.order == 2 ? @sprintf(", min scaled Jacobian %.3f", s.min_jacobian) : "",
                s.ninverted > 0 ? ", $(s.ninverted) INVERTED" : "")
    end
    if !isempty(r.interfaces)
        println(io, "  interfaces: ", join((@sprintf("%s|%s %.4g", i.phases..., i.area) for i in r.interfaces), ", "))
    end
    if !isempty(r.periodicity)
        print(io, "  periodicity:")
        for p in r.periodicity
            print(io, " ", "xyz"[p.direction], p.ok ? " ✓" : " ✗", " (", p.npairs, " pairs)")
        end
        println(io)
    end
    isempty(r.files) || print(io, "  files: ", join(r.files, ", "))
    return nothing
end
