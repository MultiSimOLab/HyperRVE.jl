"""Volume and identification of one phase in a generated RVE."""
struct PhaseInfo
    name::String
    tag::Int
    entities::Vector{Int}
    volume::Float64
    fraction::Float64
end

"""A physical group of the generated model."""
struct GroupInfo
    name::String
    dim::Int
    tag::Int
    entities::Vector{Int}
end

"""Summary statistics of the generated mesh (lengths in the units of the cell)."""
struct MeshStats
    nnodes::Int
    nelements::Int
    order::Int
    min_quality::Float64
    mean_quality::Float64
    hmin::Float64
    hmax::Float64
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
    stats::Union{Nothing,MeshStats}
    periodicity::Vector{PeriodicityCheck}
    files::Vector{String}
    model::String
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
computed from the exact geometry.
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
                             s.min_quality, s.hmin, s.hmax))
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
