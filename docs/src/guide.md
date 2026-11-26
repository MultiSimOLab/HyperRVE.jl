# User guide

```@meta
CurrentModule = RVEs
```

## Workflow

1. Describe the **cell** with [`Cell`](@ref): size, origin and periodic directions.
2. Describe the **phases** with [`Phase`](@ref): a name and one or more geometries built from
   primitives, booleans and transforms ([`Geometry`](@ref)).
3. Wrap them in an [`RVE`](@ref) (background `matrix`, overlap rule, supercell tiling).
4. Call [`generate`](@ref), which returns an [`RVEMesh`](@ref) and writes the mesh files.

```julia
using RVEs
cell = Cell((1.0, 1.0, 1.0))
rve  = RVE(cell, Phase(:particles, [Sphere((0.5, 0.5, 0.5), 0.25), Sphere((0, 0, 0), 0.25)]; refine = 0.03))
res  = generate(rve; h = 0.1, output = "bcc.msh")
volume_fraction(res, :particles)
```

## Cells

`Cell(size; origin, periodic)` accepts 2 or 3 components; lengths are unit-agnostic and go to
the mesh exactly as given (mm in, mm out). `cellpoint(cell, ξ)` maps fractional coordinates to
physical ones. Cells of any size work: the geometry is built and meshed in a normalised frame
and mapped back, so a 1 µm cell in metres gives the same mesh as a unit cell.

`RVE(...; tiling = (n1, n2, n3))` or [`supercell`](@ref) repeat the cell content.

## Geometry

| 3D | 2D |
|----|----|
| `Sphere(center, radius)` | `Disk(center, radius)` |
| `Ellipsoid(center, radii)` | `Ellipse(center, radii)` |
| `Box(corner, size)` / `Box(center = c, size = s)` | `Rectangle(corner, size)` |
| `Cylinder(base, axis, radius)` (finite) | |
| `Capsule(p1, p2, radius)` (cylinder with spherical caps) | |
| `Fibre(point, direction, radius)` (infinite) | |

Booleans: `Fuse(a, b, ...)` (`a ∪ b`), `Cut(object, tools...)` (`setdiff`), `Intersect(a, b, ...)`
(`a ∩ b`). Transforms: `Translate(g, shift)`, `Rotate(g; angle, axis, center)` (about the
geometry's own centre by default), `Scale(g, factors; center)`.

Every item of a phase that crosses a periodic face, edge or corner is completed with its
periodic images; items larger than the cell are handled too. Fibres must be parallel to a
cell axis to be phase items (oblique fibres can be used inside booleans).

## Phases and overlaps

`Phase(name, items; fuse = true)` unites the items of a phase and their periodic images into
one body before meshing (no internal surfaces where items overlap): use it for strut
lattices and clusters of overlapping particles.

Phases are painted in order: where two phases overlap, the one listed later wins
(`overlaps = :paint`); `overlaps = :error` makes any overlap an error. A coated particle is
`[Phase(:coating, Sphere(c, R)), Phase(:core, Sphere(c, r))]`. `matrix = nothing` removes the
background (porous solids, lattices): the mesh then contains only the phases.

## Mesh control

- `generate(rve; h)`: bulk size; `Phase(...; h)`: size inside a phase.
- `Phase(...; refine = 0.02)`: size on the interfaces of the phase, kept over `width` and then
  growing with ratio `growth` ([`InterfaceRefinement`](@ref), [`MeshOptions`](@ref)).
- The size field is periodic: refining an inclusion close to a face refines the opposite face.
- [`MeshOptions`](@ref): `hmin`, `hmax`, `curvature`, `order` (1 or 2, curved or straight),
  `algorithm`, `surface_algorithm`, `optimize`, `optimize_curved`.
- [`quality_report`](@ref) prints the quality histogram, the minimum per phase and the worst
  element; `generate` warns when the quality is poor or curved elements are inverted.

## Physical groups

| name | dimensions | tag |
|------|-----------|-----|
| matrix, phase names | D | 1, 1 + i |
| `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax` | D-1 … 0 | 1–6 |
| `corners`, `boundary` | 0 / D-1 … 0 | 7, 8 |
| `edge_<f1>_<f2>` (3D) | 1, 0 | 11–22 |
| `corner_x…_y…[_z…]` | 0 | 31–38 |
| `interfaces`, `interface_<a>_<b>` | D-1 … 0 | 99, from 102 |

The result object gives the same information: [`group_names`](@ref), [`group_tag`](@ref),
[`interface_names`](@ref), [`interface_area`](@ref), [`interface_tag`](@ref).

## Using the meshes in HyperFEM / Gridap

```julia
using Gridap, GridapGmsh
model = GmshDiscreteModel("rve.msh")              # periodic nodes are identified
Ωp = Triangulation(model, tags = "particles")
V  = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1); dirichlet_tags = ["corners"])
Γ  = BoundaryTriangulation(model; tags = "interface_matrix_particles")
```

- In periodic directions GridapGmsh merges the labels of opposite faces: use `corners` and the
  phase names there; the face labels are reliable in non-periodic directions.
- GridapGmsh 0.7 reads periodic links only for first-order meshes: with periodic cells use
  `order = 1` and higher-order finite element spaces.
- `BoundaryTriangulation(model; tags = "interface_a_b")` takes the cells of a side chosen by
  Gridap (normal orientation not controlled); for one-sided traces use
  `BoundaryTriangulation(Triangulation(model, tags = "a"); tags = "interface_a_b")`, and for
  jumps `InterfaceTriangulation(Triangulation(model, tags = "b"), Triangulation(model, tags = "a"))`.

## Known limitations

- Bodies flush with a periodic face: their boundary on the face is treated as cell face.
- Gmsh may fail to mesh periodic curves for bodies in some generic positions very close to a
  cell edge or corner; the error is explicit.
- Copies of items much larger than the cell make the geometry step slow.
- Volumes and areas of pieces bounded by B-spline surfaces (ellipsoids and bodies transformed
  with `Scale`) are accurate to about 1e-2.
- Gmsh's default surface mesher sometimes leaves such B-spline surfaces without elements;
  `generate` remeshes those surfaces with MeshAdapt, and stops with an error if a phase is
  still not meshed completely.
