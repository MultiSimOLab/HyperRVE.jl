# Gmsh session ownership. `generate` initialises Gmsh only if nobody did, and finalises
# whatever it initialised (always on error, unless `keep_open` on success).

_gmsh_isopen() = Bool(gmsh.isInitialized())

_nopopup() = "-nopopup" in ARGS || get(ENV, "RVES_NOPOPUP", "") == "1"

function _with_session(f, name::AbstractString; keep_open::Bool, verbosity::Int)
    owned = !_gmsh_isopen()
    owned && gmsh.initialize()
    ok = false
    try
        gmsh.option.setNumber("General.Terminal", verbosity >= 2 ? 1 : 0)
        gmsh.model.add(name)
        r = f()
        ok = true
        return r
    catch
        if !owned && _gmsh_isopen() && name in gmsh.model.list()
            gmsh.model.setCurrent(name)
            gmsh.model.remove()
        end
        rethrow()
    finally
        if owned && !(keep_open && ok) && _gmsh_isopen()
            gmsh.finalize()
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
`generate(...; keep_open = true)`). The format follows the extension (`.msh`, `.vtk`, ...).
"""
function write_mesh(r::RVEMesh, path::AbstractString)
    _require_live(r)
    p = abspath(path)
    gmsh.write(p)
    push!(r.files, p)
    return p
end

"""
    view_mesh(res)

Open the Gmsh GUI on a live result (blocking). Ignored when Julia runs with `-nopopup` in
`ARGS` or with `ENV["RVES_NOPOPUP"] = "1"`.
"""
function view_mesh(r::RVEMesh)
    _require_live(r)
    _nopopup() ? (@info "view_mesh: GUI disabled (-nopopup / RVES_NOPOPUP)") : gmsh.fltk.run()
    return nothing
end

"""
    close_gmsh()

Finalise Gmsh if it is initialised (for sessions left open with `keep_open = true`).
"""
close_gmsh() = (_gmsh_isopen() && gmsh.finalize(); nothing)
