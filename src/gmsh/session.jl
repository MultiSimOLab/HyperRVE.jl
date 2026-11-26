# Gmsh session ownership. `generate` initialises Gmsh only if nobody did, and finalises
# whatever it initialised (always on error, unless `keep_open` on success).

_gmsh_isopen() = Bool(gmsh.isInitialized())

_nopopup() = "-nopopup" in ARGS || get(ENV, "RVES_NOPOPUP", "") == "1"

# Options changed by `generate`; restored afterwards when the session belongs to the caller.
const _SESSION_OPTIONS = ("General.Terminal", "Geometry.OCCBoundsUseStl", "Mesh.MeshSizeExtendFromBoundary",
                          "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature",
                          "Mesh.MshFileVersion", "Mesh.SaveAll", "Mesh.Binary")

function _unique_model_name(name)
    names = gmsh.model.list()
    name in names || return String(name)
    k = 2
    while "$(name)_$k" in names
        k += 1
    end
    return "$(name)_$k"
end

function _with_session(f, name::AbstractString; keep_open::Bool, verbosity::Int)
    owned = !_gmsh_isopen()
    # user configuration files (~/.gmshrc, ~/.gmsh-options) are not read: results must not
    # depend on them (they can, e.g., switch the output to MSH 2.2)
    owned && gmsh.initialize(String[], false)
    saved = owned ? nothing : [(o, gmsh.option.getNumber(o)) for o in _SESSION_OPTIONS]
    previous = owned ? "" : gmsh.model.getCurrent()
    model = owned ? String(name) : _unique_model_name(name)
    ok = false
    try
        gmsh.option.setNumber("General.Terminal", verbosity >= 2 ? 1 : 0)
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.model.add(model)
        r = f(model)
        ok = true
        return r
    catch
        if !owned && _gmsh_isopen() && model in gmsh.model.list()
            gmsh.model.setCurrent(model)
            gmsh.model.remove()
        end
        rethrow()
    finally
        if owned && !(keep_open && ok) && _gmsh_isopen()
            gmsh.finalize()
        elseif saved !== nothing && _gmsh_isopen()
            for (o, v) in saved
                gmsh.option.setNumber(o, v)
            end
            previous in gmsh.model.list() && gmsh.model.setCurrent(previous)
        end
    end
end

function _require_live(r)
    if !_gmsh_isopen() || !(r.model in gmsh.model.list())
        error("the Gmsh model of this result is closed; pass output = ... to generate, " *
              "or use `generate(rve) do res ... end` / `keep_open = true`")
    end
    gmsh.model.setCurrent(r.model)
    return nothing
end

"""
    write_mesh(res, path)

Write the mesh of a live result (inside `generate(rve) do res ... end`, or after
`generate(...; keep_open = true)`). The format follows the extension (`.msh`, `.vtk`, ...);
`.msh` and `.msh4` files are written in MSH 4.1 ASCII with the periodic links in physical
units (other Gmsh extensions such as `.msh2` select other formats and are written as Gmsh
does, without that correction).

Use `write_mesh` rather than `gmsh.write`: the model is built in a normalised frame (the CAD
entities stay in it, only the mesh nodes are mapped back to physical coordinates), so a
`.msh` written directly with `gmsh.write` stores periodic maps in frame units and Gmsh would
misplace the periodic nodes when reading it back.
"""
function write_mesh(r::RVEMesh, path::AbstractString)
    previous = _gmsh_isopen() ? gmsh.model.getCurrent() : ""
    _require_live(r)
    p = abspath(path)
    isdir(dirname(p)) || throw(ArgumentError("output directory $(dirname(p)) does not exist"))
    v = gmsh.option.getNumber("Mesh.MshFileVersion")
    b = gmsh.option.getNumber("Mesh.Binary")
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    try
        gmsh.write(p)
        _physical_periodic_affine!(p, Frame(r.rve))
    finally
        gmsh.option.setNumber("Mesh.MshFileVersion", v)
        gmsh.option.setNumber("Mesh.Binary", b)
        previous in gmsh.model.list() && gmsh.model.setCurrent(previous)
    end
    p in r.files || push!(r.files, p)
    return p
end

"""
    view_mesh(res)

Open the Gmsh GUI on a live result (blocking). Ignored when Julia runs with `-nopopup` in
`ARGS` or with `ENV["RVES_NOPOPUP"] = "1"`.
"""
function view_mesh(r::RVEMesh)
    _require_live(r)
    _run_gui()
    return nothing
end

"""
    close_gmsh()

Finalise Gmsh if it is initialised (for sessions left open with `keep_open = true`).
"""
close_gmsh() = (_gmsh_isopen() && gmsh.finalize(); nothing)

# The CAD model stays in normalised coordinates while the mesh is physical: show only the mesh.
function _run_gui()
    if _nopopup()
        @info "Gmsh GUI disabled (-nopopup / RVES_NOPOPUP)"
        return nothing
    end
    for o in ("Geometry.Points", "Geometry.Curves", "Geometry.Surfaces", "Geometry.Volumes")
        gmsh.option.setNumber(o, 0)
    end
    gmsh.fltk.run()
    return nothing
end
