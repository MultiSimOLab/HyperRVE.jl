"""
    generate(rve; h = nothing, mesh = MeshOptions(), output = nothing, geometry_only = false,
             check_periodicity = true, view = false, keep_open = false, verbosity = 0,
             name = "RVE") -> RVEMesh
    generate(f, rve; kwargs...)

Build the geometry of `rve` in Gmsh, assign phases, create physical groups, impose
periodicity and mesh it. Gmsh is initialised and finalised automatically; the returned
[`RVEMesh`](@ref) is plain Julia data.

- `h`: global element size (shorthand for `mesh = MeshOptions(h = h)`).
- `output`: file name or vector of file names written before Gmsh is closed (`.msh`, `.vtk`...).
- `geometry_only`: stop before meshing (volumes and groups only; fast).
- `check_periodicity`: verify node by node that opposite periodic faces have matching meshes
  (throws [`PeriodicityError`](@ref) otherwise).
- `view`: open the Gmsh GUI before returning (skipped with `-nopopup` or `RVES_NOPOPUP=1`).
- `keep_open`: leave Gmsh initialised so that [`write_mesh`](@ref) and [`view_mesh`](@ref)
  can be used; close it with [`close_gmsh`](@ref).
- `verbosity`: 0 silent, 1 summary, 2 Gmsh terminal output.

The do-block form runs `f(res)` while the Gmsh model is alive, then closes Gmsh.

Physical groups: one per phase (dimension 3; matrix tag 1, phase `i` tag `1 + i`), and on the
cell boundary `xmin, xmax, ymin, ymax, zmin, zmax` (tags 1-6, closed faces, defined in
dimensions 2, 1 and 0), `corners` (7), `boundary` (8), `edge_<face>_<face>` (11-22) and
`corner_<x>_<y>_<z>` (31-38).
"""
generate(rve::RVE; kwargs...) = _generate(nothing, rve; kwargs...)
generate(f::Function, rve::RVE; kwargs...) = _generate(f, rve; kwargs...)

function _generate(f, rve::RVE{D}; h=nothing, mesh::MeshOptions=MeshOptions(), output=nothing,
                   geometry_only::Bool=false, check_periodicity::Bool=true, view::Bool=false,
                   keep_open::Bool=false, verbosity::Integer=0, name::AbstractString="RVE") where {D}
    if h !== nothing
        mesh.h === nothing || throw(ArgumentError("give the global size either as h or as mesh.h, not both"))
        mesh = MeshOptions(; h=_positive(h, "h"), order=mesh.order, algorithm=mesh.algorithm,
                           optimize=mesh.optimize, recombine=mesh.recombine)
    end
    _check(mesh)
    files = output === nothing ? String[] :
            output isa AbstractString ? [abspath(output)] : [abspath(String(o)) for o in output]
    frame = Frame(rve)
    return _with_session(name; keep_open=keep_open && f === nothing, verbosity=Int(verbosity)) do
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
        part = assemble!(rve, frame)
        groups = add_physical_groups!(rve, frame, part)
        set_periodicity!(rve, frame)
        vdom = volume(rve)
        phases = PhaseInfo[]
        for (n, t, e) in zip(part.names, part.tags, part.entities)
            v = sum((from_frame_measure(frame, gmsh.model.occ.getMass(D, x)) for x in e); init=0.0)
            push!(phases, PhaseInfo(n, isempty(e) ? 0 : t, Int.(e), v, v / vdom))
        end
        stats = nothing
        checks = PeriodicityCheck[]
        if !geometry_only
            mesh!(rve, frame, part, mesh)
            restore_frame!(frame)
            check_periodicity && (checks = verify_periodicity(rve, frame))
            stats = mesh_stats(frame, mesh)
        end
        for p in files
            gmsh.write(p)
        end
        res = RVEMesh{D}(rve, mesh, phases, groups, stats, checks, files, String(name))
        verbosity >= 1 && @info sprint(show, MIME("text/plain"), res)
        f === nothing || f(res)
        if view
            _nopopup() ? (@info "generate: GUI disabled (-nopopup / RVES_NOPOPUP)") : gmsh.fltk.run()
        end
        res
    end
end
