"""
    generate(rve; h = nothing, mesh = MeshOptions(), output = nothing, geometry_only = false,
             check_periodicity = true, view = false, keep_open = false, verbosity = 0,
             name = "RVE") -> RVEMesh
    generate(f, rve; kwargs...)

Build the geometry of `rve` in Gmsh, assign phases, create physical groups, impose
periodicity and mesh it. Gmsh is initialised and finalised automatically; the returned
[`RVEMesh`](@ref) is plain Julia data.

- `h`: global element size (shorthand for `mesh = MeshOptions(h = h)`); it is reduced near
  inclusions with refinement and near boundaries of small curvature radius (see [`MeshOptions`](@ref)).
- `output`: file name or vector of file names written before Gmsh is closed (`.msh`, `.vtk`...).
- `geometry_only`: stop before meshing (volumes and groups only; fast).
- `check_periodicity`: verify node by node that opposite periodic faces have matching meshes
  (throws [`PeriodicityError`](@ref) otherwise).
- `view`: open the Gmsh GUI before returning (skipped with `-nopopup` or `RVES_NOPOPUP=1`).
- `keep_open`: leave Gmsh initialised so that [`write_mesh`](@ref) and [`view_mesh`](@ref)
  can be used; close it with [`close_gmsh`](@ref).
- `verbosity`: 0 silent, 1 summary, 2 Gmsh terminal output.
- `name`: name of the Gmsh model (made unique if a model with that name is already open).
- `interfaces`: create the interface physical groups (see below).

The do-block form runs `f(res)` while the Gmsh model is alive, then closes Gmsh (unless
`keep_open = true`); it returns `res`, not the value of `f`. In the live model the mesh is in
physical coordinates but the CAD entities are in a normalised frame (domain origin at 0,
longest edge 1): use [`write_mesh`](@ref), not `gmsh.write`, to export it. Each call adds a
new Gmsh model to an open session (names made unique); call [`close_gmsh`](@ref) to free them.

Gmsh is initialised without reading the user configuration files, and `.msh` files are
written in format 4.1. If the caller already initialised Gmsh, that session is left open and
the options changed by `generate` are restored on return.

Physical groups: one per phase (dimension `D`; matrix tag 1, phase `i` tag `1 + i`), and on
the cell boundary `xmin, xmax, ymin, ymax, zmin, zmax` (tags 1-6, closed faces, defined in
every lower dimension), `corners` (7), `boundary` (8), `edge_<face>_<face>` (11-22, 3D only)
and `corner_<x>_<y>[_<z>]` (31-38). Interfaces between phases (surfaces not on the cell faces
whose two sides belong to different phases; one side `void` for free surfaces without
matrix) get the groups `interface_<a>_<b>` (phases in the order of `phase_names(res)`; with
`N` phases including the matrix and `i < j` their positions, `void` being `N + 1`, the tag is
`100 + (i-1)(N+1) + j`) and their union `interfaces` (tag 99), defined in dimensions D-1 to 0
(closed sets; closed surfaces such as a whole sphere contribute their seam curve and poles).
In Gridap, `BoundaryTriangulation(model; tags = "interface_a_b")` gives the interface
facets with the cells of one side chosen by Gridap (the orientation of the normal is not
controlled: fine for areas, Dirichlet data and integrals of scalar facet quantities). For a
one-sided trace use `BoundaryTriangulation(Triangulation(model, tags = "a"); tags =
"interface_a_b")`; for jumps and normals oriented from one phase to the other use
`InterfaceTriangulation(Triangulation(model, tags = "b"), Triangulation(model, tags = "a"))`
(plus side in `b`), which does not need the interface groups. Facets on an interface carry
the interface labels instead of the phase labels. In 2D the faces are the four sides of the cell (curves
and points) and the mesh is made of triangles in the plane z = 0. Groups without entities are
not created: with `matrix = nothing` the cell faces exist only where the phases touch them, and
a periodic direction without any entity on its faces is reported with 0 pairs.

When a periodic mesh is read by GridapGmsh, master and slave entities are identified, so the
labels `xmin`/`xmax` (and the `corner_*`, `edge_*` groups) of a periodic direction are shared
between both sides and are not reliable as boundary tags; in periodic directions use
`corners` (e.g. to fix rigid-body motions) and the phase names.
"""
generate(rve::RVE; kwargs...) = _generate(nothing, rve; kwargs...)
generate(f::Function, rve::RVE; kwargs...) = _generate(f, rve; kwargs...)

function _generate(f, rve::RVE{D}; h=nothing, mesh::MeshOptions=MeshOptions(), output=nothing,
                   geometry_only::Bool=false, check_periodicity::Bool=true, view::Bool=false,
                   keep_open::Bool=false, verbosity::Integer=0, name::AbstractString="RVE",
                   interfaces::Bool=true) where {D}
    if h !== nothing
        mesh.h === nothing || throw(ArgumentError("give the global size either as h or as mesh.h, not both"))
        mesh = _with(mesh; h=_positive(h, "h"))
    end
    mesh.optimize_curved == :optimize && any(rve.cell.periodic) &&
        throw(ArgumentError("optimize_curved = :optimize moves nodes of periodic faces; use :elastic for periodic cells"))
    files = output === nothing ? String[] :
            output isa AbstractString ? [abspath(output)] : [abspath(String(o)) for o in output]
    for p in files
        isdir(dirname(p)) || throw(ArgumentError("output directory $(dirname(p)) does not exist"))
    end
    geometry_only && !isempty(files) &&
        @warn "generate: geometry_only = true, the files written contain no mesh"
    frame = Frame(rve)
    return _with_session(name; keep_open=keep_open, verbosity=Int(verbosity)) do model
        part = assemble!(rve, frame)
        on = _boundary_entities(rve, frame)
        groups = add_physical_groups!(rve, frame, part, on)
        set_periodicity!(rve, frame, on)
        ifs = classify_interfaces(part, on)
        iinfo = interface_groups!(groups, rve, frame, part, ifs; create=interfaces)
        vdom = volume(rve)
        phases = PhaseInfo[]
        for (n, t, e) in zip(part.names, part.tags, part.entities)
            v = sum((from_frame_measure(frame, gmsh.model.occ.getMass(D, x)) for x in e); init=0.0)
            push!(phases, PhaseInfo(n, isempty(e) ? 0 : t, Int.(e), v, v / vdom))
        end
        stats = nothing
        checks = PeriodicityCheck[]
        if !geometry_only
            mesh!(rve, frame, part, ifs, mesh, on)
            restore_frame!(frame)
            check_periodicity && (checks = verify_periodicity(rve, frame))
            stats = mesh_stats(frame, mesh, part)
            stats.nelements > 0 ||
                error("Gmsh produced no $(D == 3 ? "tetrahedra" : "triangles"); the mesh size h may be " *
                      "too large for the smallest inclusion (try a smaller h or InterfaceRefinement)")
            if stats.min_quality < 0.02 || stats.ninverted > 0
                @warn "generate: poor element quality (minimum SICN $(round(stats.min_quality; sigdigits=2)), " *
                      "$(stats.ninverted) inverted) in phase $(stats.worst_phase) near " *
                      "$(round.(stats.worst_point; sigdigits=4)). " *
                      (stats.order == 2 && stats.ninverted > 0 ?
                       "Curved order-2 elements are inverted where the mesh is coarse compared with the curvature: " *
                       "refine that phase, raise MeshOptions(curvature = ...), use optimize_curved = :elastic or curved = false." :
                       "Common causes: a thin ligament between inclusions or between an inclusion and a face that is " *
                       "coarser than the local size; refine that phase (refine = ..., or Phase h), increase the gap, " *
                       "reduce h, or use optimize = :netgen.")
            end
        end
        for p in files
            gmsh.write(p)
            _physical_periodic_affine!(p, frame)
        end
        res = RVEMesh{D}(rve, mesh, phases, groups, iinfo, stats, checks, files, model)
        verbosity >= 1 && @info sprint(show, MIME("text/plain"), res)
        f === nothing || f(res)
        if view
            _run_gui()
        end
        res
    end
end
