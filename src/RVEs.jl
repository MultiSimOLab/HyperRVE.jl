"""
    RVEs

Periodic representative volume element (RVE) meshes with the Gmsh Julia API.

Describe a cell ([`Cell`](@ref)), the phases inside it ([`Phase`](@ref)) built from geometric
primitives, booleans and transforms, wrap everything in an [`RVE`](@ref) and call
[`generate`](@ref). Periodic images of inclusions crossing the cell boundary, phase
assignment, physical groups and periodic meshing are handled automatically.
"""
module RVEs

using Gmsh: Gmsh, gmsh
using Printf: @sprintf
using Random: Random, AbstractRNG

export Cell, cellpoint
export Geometry, Sphere, Ellipsoid, Box, Cylinder, Capsule, Fibre, Disk, Ellipse, Rectangle, Polyhedron, Polygon
export Fuse, Cut, Intersect, Translate, Rotate, Scale
export Phase, RVE, supercell
export struts, lattice, lattice_topology, relative_density, strut_radius
export random_spheres, random_disks, random_fibres
export MeshOptions, InterfaceRefinement
export generate, RVEMesh, write_mesh, view_mesh, close_gmsh, PeriodicityError
export volume, volume_fraction, volume_fractions, phase_names, group_names, group_tag, quality_report
export interface_names, interface_area, interface_tag

include("vectors.jl")
include("cell.jl")
include("geometry.jl")
include("phase.jl")
include("meshoptions.jl")
include("rve.jl")
include("lattices.jl")
include("random.jl")
include("frame.jl")
include("result.jl")
include("gmsh/session.jl")
include("gmsh/build.jl")
include("gmsh/assemble.jl")
include("gmsh/groups.jl")
include("gmsh/periodicity.jl")
include("gmsh/interfaces.jl")
include("gmsh/meshing.jl")
include("generate.jl")

end
