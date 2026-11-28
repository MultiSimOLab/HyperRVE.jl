"""
    HyperRVE

Periodic representative volume element (RVE) meshes with the Gmsh Julia API.

Describe a cell ([`Cell`](@ref)), the phases inside it ([`Phase`](@ref)) built from geometric
primitives, booleans and transforms, wrap everything in an [`RVE`](@ref) and call
[`generate`](@ref). Periodic images of inclusions crossing the cell boundary, phase
assignment, physical groups and periodic meshing are handled automatically.
"""
module HyperRVE

using Gmsh: Gmsh, gmsh
using LinearAlgebra: LinearAlgebra
using Printf: @sprintf
using Random: Random, AbstractRNG

export Cell, cellpoint
export Geometry, Sphere, Ellipsoid, Box, Cylinder, Capsule, Fibre, Disk, Ellipse, Rectangle, Polyhedron, Polygon
export Fuse, Cut, Intersect, Translate, Rotate, Scale
export Phase, RVE, supercell
export struts, lattice, lattice_topology, relative_density, strut_radius, lattice_cell
export hollow_struts, hollow_lattice
export random_spheres, random_disks, random_fibres, random_capsules, orientation_tensor
export Tessellation, voronoi, grains, grain_phases, grain_volumes
export dense_spheres, laguerre, grain_diameters
export truncated_lognormal, TruncatedLogNormal, laguerre_foam, Foam, pore_statistics, open_pores, perforate, foam_phases
export tessellation_edges, foam_struts, pores, wall_thickness, pore_rounding
export particle_chains
export particle_array, laminate
export plate_cells, plate_density, plate_thickness
export tpms, tpms_level, tpms_density, TPMS
export chiral, rotating_squares, rotating_triangles
export MeshOptions, InterfaceRefinement, GapRefinement
export generate, RVEMesh, write_mesh, view_mesh, close_gmsh, PeriodicityError
export volume, volume_fraction, volume_fractions, phase_names, group_names, group_tag, quality_report
export interface_names, interface_area, interface_tag

include("vectors.jl")
include("cell.jl")
include("geometry.jl")
include("phase.jl")
include("gaps_types.jl")
include("meshoptions.jl")
include("rve.jl")
include("lattices.jl")
include("random.jl")
include("capsules.jl")
include("voronoi.jl")
include("foams.jl")
include("laguerre.jl")
include("rounding.jl")
include("metalfoam.jl")
include("windows.jl")
include("perforations.jl")
include("plates.jl")
include("tpms.jl")
include("chains.jl")
include("arrays.jl")
include("auxetic.jl")
include("frame.jl")
include("result.jl")
include("gmsh/session.jl")
include("gmsh/build.jl")
include("gmsh/assemble.jl")
include("gmsh/groups.jl")
include("gmsh/periodicity.jl")
include("gmsh/interfaces.jl")
include("gmsh/meshing.jl")
include("gaps.jl")
include("gmsh/tpms.jl")
include("generate.jl")

end
