# Meshing path of TPMS phases (experimental): the triangulated surface of the kernel is added as
# discrete Gmsh entities (corner points, cap patches on the faces, interface pieces), closed into
# volumes and meshed without remeshing the surface. OpenCASCADE is not involved.

_tpms_items(rve) = [it for p in rve.phases for it in p.items]

_has_tpms(::Geometry) = false
_has_tpms(::TPMS) = true
_has_tpms(g::Union{Translate,Rotate,Scale}) = _has_tpms(g.child)
_has_tpms(g::Union{Fuse,Intersect}) = any(_has_tpms, g.parts)
_has_tpms(g::Cut) = _has_tpms(g.object) || any(_has_tpms, g.tools)

function _is_tpms_rve(rve::RVE)
    items = _tpms_items(rve)
    any(_has_tpms, items) || return false
    msg = "TPMS phases are experimental and cannot be combined with other features: "
    (length(rve.phases) == 1 && length(items) == 1 && items[1] isa TPMS) ||
        throw(ArgumentError(msg * "the TPMS must be the only item of the only phase (a matrix or matrix = nothing)"))
    p = rve.phases[1]
    (p.h === nothing && p.refine === nothing && !p.fuse && !p.void) ||
        throw(ArgumentError(msg * "no h, refine, fuse or void in its phase"))
    c = rve.cell
    all(c.periodic) || throw(ArgumentError(msg * "the cell must be periodic in the three directions"))
    maximum(c.size) - minimum(c.size) <= 1e-12 * maximum(c.size) || throw(ArgumentError(msg * "the cell must be cubic"))
    all(==(1), rve.tiling) || throw(ArgumentError(msg * "no tiling; use the periods of tpms"))
    return true
end

# connected patches (by shared edges) of the triangles with label lab
function _tp_patches(S, idx)
    parent = collect(eachindex(idx))
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    owner = Dict{NTuple{2,Int},Int}()
    for (k, i) in enumerate(idx), m in 1:3
        t = S.tri[i]
        e = minmax(t[m], t[mod1(m + 1, 3)])
        j = get!(owner, e, k)
        a, b = find(k), find(j)
        a == b || (parent[a] = b)
    end
    roots = Dict{Int,Int}()
    return [get!(roots, find(k), length(roots) + 1) for k in eachindex(idx)]
end

function _tp_inside(S, p, tris)          # ray casting along a generic direction
    d = (0.5773, 0.7071, 0.4082)
    n = 0
    for i in tris
        a, b, c = S.P[S.tri[i][1]], S.P[S.tri[i][2]], S.P[S.tri[i][3]]
        e1 = b .- a; e2 = c .- a
        h = (d[2] * e2[3] - d[3] * e2[2], d[3] * e2[1] - d[1] * e2[3], d[1] * e2[2] - d[2] * e2[1])
        det = sum(e1 .* h)
        abs(det) < 1e-14 && continue
        s = p .- a
        u = sum(s .* h) / det
        (u < 0 || u > 1) && continue
        q = (s[2] * e1[3] - s[3] * e1[2], s[3] * e1[1] - s[1] * e1[3], s[1] * e1[2] - s[2] * e1[1])
        v = sum(d .* q) / det
        (v < 0 || u + v > 1) && continue
        sum(e2 .* q) / det > 1e-12 && (n += 1)
    end
    return isodd(n)
end

# builds the discrete model in the unit cube; returns volumes per phase slot (1 solid, 2 complement),
# surfaces of each cap patch and interface piece, and the corner points
function _tp_build!(S, matrix::Bool)
    model = gmsh.model
    tv(i) = (a = S.P[S.tri[i][1]]; b = S.P[S.tri[i][2]]; c = S.P[S.tri[i][3]];
             (a[1] * (b[2] * c[3] - b[3] * c[2]) + a[2] * (b[3] * c[1] - b[1] * c[3]) + a[3] * (b[1] * c[2] - b[2] * c[1])) / 6)
    cs, ns = _tp_shells(S, 0:6)
    cm, nm = matrix ? _tp_shells(S, (0, 7:12...)) : (zeros(Int, length(S.tri)), 0)
    # entities: caps split into connected patches in their plane; interface pieces by (solid shell, matrix shell)
    groups = Dict{NTuple{3,Int},Vector{Int}}()           # (label, patch or solid shell, matrix shell) -> triangles
    for lab in 1:12
        (!matrix && lab > 6) && continue
        idx = findall(==(lab), S.lab)
        isempty(idx) && continue
        for (i, p) in zip(idx, _tp_patches(S, idx))
            push!(get!(groups, (lab, p, 0), Int[]), i)
        end
    end
    for (i, l) in enumerate(S.lab)
        l == 0 && push!(get!(groups, (0, cs[i], cm[i]), Int[]), i)
    end
    owned = falses(length(S.P))
    added = Tuple{Int,Int,Vector{Int}}[]                   # (dim, entity, kernel nodes in the order added)
    # corner points first (groups corners, corner_*), then caps (they own the curves on the faces), then the interface
    # without matrix, only the solid side is meshed
    used = falses(length(S.P))
    keep(l) = matrix || l <= 6
    for (t, l) in zip(S.tri, S.lab)
        keep(l) && (used[collect(t)] .= true)
    end
    onb(x) = abs(x) < 1e-12 || abs(x - 1) < 1e-12
    oniface = falses(length(S.P))
    for (t, l) in zip(S.tri, S.lab)
        l == 0 && (oniface[collect(t)] .= true)
    end
    # points: the cell corners, and the nodes where the interface meets a cell edge (so that the
    # interface groups are closed sets, as in the OCC path)
    corners = Dict{NTuple{3,Int},Int}()
    points = Dict{Int,Int}()                               # kernel node -> point entity
    for (p, x) in enumerate(S.P)
        used[p] || continue
        nb = count(onb, x)
        (nb == 3 || (nb == 2 && oniface[p])) || continue
        t = model.addDiscreteEntity(0)
        gmsh.model.mesh.addNodes(0, t, [p], collect(x))
        push!(added, (0, t, [p]))
        owned[p] = true
        points[p] = t
        nb == 3 && (corners[round.(Int, x)] = t)
    end
    # cell edges: (axis a, fixed values of the other two coordinates) -> discrete curve with the
    # triangle edges that lie on it
    segs = Dict{Tuple{Int,NTuple{2,Int}},Set{NTuple{2,Int}}}()
    for (t, l) in zip(S.tri, S.lab), m in 1:3
        keep(l) || continue
        a, b = t[m], t[mod1(m + 1, 3)]
        xa, xb = S.P[a], S.P[b]
        # both ends on the same cell edge (corners belong to every edge through them)
        for ax in 1:3
            o = [i for i in 1:3 if i != ax]
            all(i -> onb(xa[i]) && onb(xb[i]) && abs(xa[i] - xb[i]) < 1e-12, o) || continue
            push!(get!(segs, (ax, Tuple(round(Int, xa[i]) for i in o)), Set{NTuple{2,Int}}()), minmax(a, b))
        end
    end
    edges = Dict{Tuple{Int,NTuple{2,Int}},Int}()
    for k in sort!(collect(keys(segs)))
        t = model.addDiscreteEntity(1)
        edges[k] = t
        sg = sort!(collect(segs[k]))
        nodes = unique(Iterators.flatten(sg))
        new = [p for p in nodes if !owned[p]]
        owned[new] .= true
        isempty(new) || gmsh.model.mesh.addNodes(1, t, new, collect(Iterators.flatten(S.P[p] for p in new)))
        push!(added, (1, t, new))
        gmsh.model.mesh.addElementsByType(t, 1, Int[], collect(Iterators.flatten(sg)))
    end
    # interface ∩ cell faces: one discrete curve per face with the interface triangle edges on it
    fsegs = Dict{Int,Set{NTuple{2,Int}}}()
    for (t, l) in zip(S.tri, S.lab), m in 1:3
        l == 0 || continue
        a, b = t[m], t[mod1(m + 1, 3)]
        for d in 1:3, side in 0:1
            abs(S.P[a][d] - side) < 1e-12 && abs(S.P[b][d] - side) < 1e-12 &&
                push!(get!(fsegs, 2d - 1 + side, Set{NTuple{2,Int}}()), minmax(a, b))
        end
    end
    fcurves = Dict{Int,Int}()
    for f in sort!(collect(keys(fsegs)))
        t = model.addDiscreteEntity(1)
        fcurves[f] = t
        sg = sort!(collect(fsegs[f]))
        new = [p for p in unique(Iterators.flatten(sg)) if !owned[p]]
        owned[new] .= true
        isempty(new) || gmsh.model.mesh.addNodes(1, t, new, collect(Iterators.flatten(S.P[p] for p in new)))
        push!(added, (1, t, new))
        gmsh.model.mesh.addElementsByType(t, 1, Int[], collect(Iterators.flatten(sg)))
    end
    surf = Dict{NTuple{3,Int},Int}()
    for k in sort!(collect(keys(groups)); by=k -> (k[1] == 0 ? 99 : k[1], k[2], k[3]))
        s = model.addDiscreteEntity(2)
        surf[k] = s
        tri = groups[k]
        nodes = unique(Iterators.flatten(S.tri[i] for i in tri))
        new = [p for p in nodes if !owned[p]]
        owned[new] .= true
        isempty(new) || gmsh.model.mesh.addNodes(2, s, new, collect(Iterators.flatten(S.P[p] for p in new)))
        push!(added, (2, s, new))
        gmsh.model.mesh.addElementsByType(s, 2, Int[], collect(Iterators.flatten(S.tri[i] for i in tri)))
    end
    shell_of(k, slot) = k[1] == 0 ? (slot == 1 ? k[2] : k[3]) :
                        (slot == 1 && k[1] <= 6) ? cs[groups[k][1]] : (slot == 2 && k[1] > 6) ? cm[groups[k][1]] : 0
    vols = (Int[], Int[])
    # shells of negative volume are cavities (islands of the other phase): inner loops of the shell around them
    for (slot, comp, nc, sgn) in ((1, cs, ns, 1), (2, cm, nm, -1))
        tris = [findall(==(c), comp) for c in 1:nc]
        vc = [sum((S.lab[i] == 0 ? sgn : 1) * tv(i) for i in tris[c]) for c in 1:nc]
        outer = [c for c in 1:nc if vc[c] > 0]
        loopof(c) = model.geo.addSurfaceLoop([s for (k, s) in surf if shell_of(k, slot) == c])
        loops = Dict(c => [loopof(c)] for c in outer)
        for c in 1:nc
            vc[c] > 0 && continue
            p = S.P[S.tri[tris[c][1]][1]]
            hosts = length(outer) == 1 ? outer : [o for o in outer if _tp_inside(S, p, tris[o])]
            isempty(hosts) && error("TPMS: a cavity of the surface lies in no closed shell; try another resolution or level")
            host = first(hosts)
            push!(loops[host], loopof(c))
        end
        for c in outer
            push!(vols[slot], model.geo.addVolume(loops[c]))
        end
    end
    model.geo.synchronize()
    return vols, surf, corners, groups, added, edges, fcurves, points
end

# $Periodic section, as Gmsh writes it: every node is a slave once, in the link of the entity that
# owns it (corner points, cell-edge curves, cap patches), paired with its translate; nodes on
# several max faces are chained through the lower directions
function _tp_periodic_links(S, surf, points, edges, fcurves, added)
    q(x) = round.(Int, x .* 1e9)
    at = Dict(q(x) => p for (p, x) in enumerate(S.P))
    owner = Dict{Int,Tuple{Int,Int}}()
    for (dim, e, ks) in added, k in ks
        owner[k] = (dim, e)
    end
    links = Tuple{Int,Int,Int,Int,Vector{Tuple{Int,Int}}}[]   # (dim, direction, slave, master, node pairs)
    # direction along which a slave entity is mapped: the first coordinate at 1 among its fixed ones
    function link!(dim, ent, d, nodes)
        e = ntuple(i -> i == d ? 1.0 : 0.0, 3)
        pairs = Tuple{Int,Int}[]
        ms = Int[]
        for p in nodes
            m = get(at, q(S.P[p] .- e), nothing)
            m === nothing && throw(PeriodicityError(d, "TPMS: node $(S.P[p]) has no translate on the opposite face"))
            push!(pairs, (p, m))
            push!(ms, owner[m][2])
        end
        isempty(pairs) && return
        push!(links, (dim, d, ent, argmax(x -> count(==(x), ms), unique(ms)), pairs))
    end
    for (p, t) in sort!(collect(points))
        d = findfirst(x -> abs(x - 1) < 1e-12, S.P[p])
        d === nothing || link!(0, t, d, [p])
    end
    for ((ax, v), t) in sort!(collect(edges))
        o = [i for i in 1:3 if i != ax]
        k = findfirst(==(1), v)
        k === nothing && continue
        link!(1, t, o[k], [p for (dim, e, ks) in added if dim == 1 && e == t for p in ks])
    end
    for (f, t) in sort!(collect(fcurves))
        iseven(f) && link!(1, t, f ÷ 2, [p for (dim, e, ks) in added if dim == 1 && e == t for p in ks])
    end
    for (k, t) in sort!(collect(surf))
        k[1] == 0 && continue
        f = mod1(k[1], 6)
        iseven(f) || continue
        link!(2, t, f ÷ 2, [p for (dim, e, ks) in added if dim == 2 && e == t for p in ks])
    end
    return links
end

# periodic links of the live TPMS models, for write_mesh
const _TP_LINKS = Dict{String,Any}()

function _tp_write_periodic!(path, links)
    any(e -> endswith(lowercase(path), e), (".msh", ".msh4")) || return nothing
    lines = readlines(path)
    i = findfirst(==("\$EndElements"), lines)
    i === nothing && return nothing
    j = findfirst(==("\$Periodic"), lines)
    if j !== nothing                                  # replace a section written by Gmsh
        k = findfirst(==("\$EndPeriodic"), lines)
        deleteat!(lines, j:k)
    end
    sec = ["\$Periodic", string(length(links))]
    for (dim, d, s, m, pairs) in links
        T = _translation(d, 1.0)
        push!(sec, "$dim $s $m")
        push!(sec, "16 " * join((repr(x) for x in T), " "))
        push!(sec, string(length(pairs)))
        append!(sec, ["$a $b" for (a, b) in pairs])
    end
    push!(sec, "\$EndPeriodic")
    i = findfirst(==("\$EndElements"), lines)
    splice!(lines, i+1:i, sec)
    open(path, "w") do io
        foreach(l -> println(io, l), lines)
    end
    return nothing
end

function _generate_tpms(f, rve::RVE{3}; h=nothing, mesh::MeshOptions=MeshOptions(), output=nothing,
                        geometry_only::Bool=false, check_periodicity::Bool=true, view::Bool=false,
                        keep_open::Bool=false, verbosity::Integer=0, name::AbstractString="RVE", interfaces=true)
    (interfaces === true || interfaces === false || interfaces === :union) || throw(ArgumentError("interfaces must be true, false or :union, got $(repr(interfaces))"))
    mesh.order == 1 || throw(ArgumentError("TPMS phases are experimental: order 1 only"))
    if h !== nothing
        mesh.h === nothing || throw(ArgumentError("give the global size either as h or as mesh.h, not both"))
        mesh = _with(mesh; h=_positive(h, "h"))
    end
    files = output === nothing ? String[] : output isa AbstractString ? [abspath(output)] : [abspath(String(o)) for o in output]
    for p in files
        isdir(dirname(p)) || throw(ArgumentError("output directory $(dirname(p)) does not exist"))
        any(e -> endswith(lowercase(p), e), (".msh", ".msh4")) ||
            @warn "generate: the periodic node pairs of a TPMS mesh are written only to .msh files (MSH 4.1); $p has none"
    end
    geometry_only && !isempty(files) && (@warn "generate: geometry_only = true, no file is written for a TPMS phase"; files = String[])
    d0 = MeshOptions()
    ignored = [f for f in fieldnames(MeshOptions) if !(f in (:h, :order)) && getfield(mesh, f) != getfield(d0, f)]
    isempty(ignored) || @warn "generate: a TPMS mesh uses only h of the MeshOptions; ignored: $(join(ignored, ", "))"
    g = only(_tpms_items(rve))
    N = g.periods * g.resolution
    S = _tpms_checked_surface(g.kind, g.type, g.level, g.periods, N, "generate")
    frame = Frame(rve)
    L = frame.scale
    matrix = rve.matrix !== nothing
    pname = rve.phases[1].name
    names = matrix ? [rve.matrix, pname] : [pname]
    vs = _tp_solid_volume(S)
    phases = PhaseInfo[]
    matrix && push!(phases, PhaseInfo(rve.matrix, 1, Int[], (1 - vs) * L^3, 1 - vs))
    push!(phases, PhaseInfo(pname, 2, Int[], vs * L^3, vs))
    iname = "interface_" * (matrix ? rve.matrix * "_" * pname : pname * "_void")
    iarea = _tp_area(S) * L^2
    return _with_session(name; keep_open=keep_open, verbosity=Int(verbosity)) do model
        groups = GroupInfo[]
        infos = InterfaceInfo[]
        stats = nothing
        checks = PeriodicityCheck[]
        if geometry_only
            push!(infos, InterfaceInfo(iname, matrix ? (rve.matrix, pname) : (pname, "void"), 0, Int[], iarea))
        else
            vols, surf, corners, tgroups, added, cedges, fcurves, points = _tp_build!(S, matrix)
            ents = matrix ? [Int32.(vols[2]), Int32.(vols[1])] : [Int32.(vols[1])]
            tags = matrix ? [1, 2] : [2]
            for (n, t, e) in zip(names, tags, ents)
                _add_group!(groups, 3, e, t, n)
            end
            phases = [PhaseInfo(p.name, p.tag, Int.(e), p.volume, p.fraction) for (p, e) in zip(phases, ents)]
            caps(sel) = [s for (k, s) in surf if k[1] != 0 && sel(mod1(k[1], 6))]
            for d in 1:3, side in 1:2
                fc = _face_tag(d, side)
                _add_group!(groups, 2, caps(==(fc)), fc, _face_name(d, side))
                _add_group!(groups, 1, [[t for ((ax, v), t) in cedges if ax != d && v[findfirst(==(d), [i for i in 1:3 if i != ax])] == side - 1];
                                        [t for (f, t) in fcurves if f == fc]], fc, _face_name(d, side))
                _add_group!(groups, 0, [t for (p, t) in points if abs(S.P[p][d] - (side - 1)) < 1e-12], fc, _face_name(d, side))
            end
            _add_group!(groups, 0, collect(values(corners)), 7, "corners")
            _add_group!(groups, 2, caps(_ -> true), 8, "boundary")
            _add_group!(groups, 1, [collect(values(cedges)); collect(values(fcurves))], 8, "boundary")
            _add_group!(groups, 0, collect(values(points)), 8, "boundary")
            etag = 11
            for a in 1:3
                d1, d2 = Tuple(setdiff(1:3, a))
                for s1 in 1:2, s2 in 1:2
                    k = (a, (s1 - 1, s2 - 1))
                    nm = "edge_" * _face_name(d1, s1) * "_" * _face_name(d2, s2)
                    haskey(cedges, k) && _add_group!(groups, 1, [cedges[k]], etag, nm)
                    _add_group!(groups, 0, [t for (p, t) in points if abs(S.P[p][d1] - (s1 - 1)) < 1e-12 && abs(S.P[p][d2] - (s2 - 1)) < 1e-12], etag, nm)
                    etag += 1
                end
            end
            tag = 31
            for idx in Iterators.product(1:2, 1:2, 1:2)
                c = idx .- 1
                haskey(corners, c) && _add_group!(groups, 0, [corners[c]], tag,
                                                 "corner_" * join((_face_name(d, idx[d]) for d in 1:3), "_"))
                tag += 1
            end
            isurf = [s for (k, s) in surf if k[1] == 0]
            ifn = Set(p for (t, l) in zip(S.tri, S.lab) if l == 0 for p in t)
            ipts = [t for (p, t) in points if p in ifn]
            for (dim, ents) in ((2, isurf), (1, collect(values(fcurves))), (0, ipts))
                interfaces === true && _add_group!(groups, dim, ents, 102, iname)
                interfaces !== false && _add_group!(groups, dim, ents, 99, "interfaces")
            end
            push!(infos, InterfaceInfo(iname, matrix ? (rve.matrix, pname) : (pname, "void"), interfaces === true ? 102 : 0, isurf, iarea))
            gmsh.option.setNumber("Mesh.MeshSizeMax", mesh.h === nothing ? 1 / N : to_frame_length(frame, mesh.h))
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
            gmsh.option.setNumber("Mesh.Optimize", 1)
            try
                gmsh.model.mesh.generate(3)
            catch e
                error("TPMS: Gmsh could not fill the volumes of the surface ($(sprint(showerror, e))); try another resolution or level")
            end
            # Gmsh renumbers the nodes when meshing: the nodes of each entity keep their order
            links = _tp_periodic_links(S, surf, points, cedges, fcurves, added)
            tagof = Dict{Int,Int}()
            for (dim, e, ks) in added
                nt = gmsh.model.mesh.getNodes(dim, e, false, false)[1]
                length(nt) == length(ks) || error("TPMS: Gmsh changed the surface nodes of entity $e")
                for (k, t) in zip(ks, nt)
                    x = gmsh.model.mesh.getNode(t)[1]
                    maximum(abs.(x .- S.P[k])) < 1e-9 || error("TPMS: Gmsh reordered the surface nodes of entity $e")
                    tagof[k] = Int(t)
                end
            end
            kt(p) = tagof[p]
            links = [(dim, d, sl, ms, [(kt(a), kt(b)) for (a, b) in pr]) for (dim, d, sl, ms, pr) in links]
            restore_frame!(frame)
            part = Partition(3, names, tags, ents, Dict{Int32,Tuple{Int,Int}}())
            stats = mesh_stats(frame, mesh, part)
            stats.min_quality < 0.02 &&
                @warn "generate: poor element quality (minimum SICN $(round(stats.min_quality; sigdigits=2))) in phase " *
                      "$(stats.worst_phase) near $(round.(stats.worst_point; sigdigits=4)); for a TPMS phase raise the resolution of tpms"
            check_completeness(frame, part, phases, mesh)
            # the mesh is the kernel's solid: its volumes must match
            for (n, e, ph) in zip(names, ents, phases)
                v = 0.0
                for x in e, t in gmsh.model.mesh.getElements(3, x)[2]
                    v += sum(gmsh.model.mesh.getElementQualities(t, "volume"))
                end
                abs(v - ph.volume) <= 1e-8 * L^3 || error("TPMS: the mesh of $n has volume $v instead of $(ph.volume); " *
                                                          "the surface could not be meshed as it is (try another resolution or level)")
            end
            check_periodicity && (checks = verify_periodicity(rve, frame))
            _TP_LINKS[model] = links
            for p in files
                gmsh.write(p)
                _tp_write_periodic!(p, links)
                _physical_periodic_affine!(p, frame)
            end
        end
        res = RVEMesh{3}(rve, mesh, phases, groups, infos, stats, checks, files, model)
        verbosity >= 1 && @info sprint(show, MIME("text/plain"), res)
        f === nothing || f(res)
        view && _run_gui()
        keep_open || delete!(_TP_LINKS, model)
        res
    end
end
