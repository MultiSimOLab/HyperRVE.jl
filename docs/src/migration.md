# Migrating from the 2024 API

version 0.2 replaces the whole public API. The table maps every old call to its replacement.

| 2024 API | 0.2 |
|----|----|
| `startGmsh()`, `stopGmsh()`, `ShowInfo(n)` | handled by `generate`; `verbosity = n` |
| `RVE(size, flags, origin, h)` | `Cell(size; origin, periodic)`, `RVE(cell, phases)`, `generate(rve; h)` |
| `Sphere(c, r, ref)` | `Sphere(c, r)` inside `Phase(name, ...; refine = InterfaceRefinement(level = ref))`, or simply `refine = <interface size>` |
| `Box(size, origin, ref)` | `Box(corner = origin, size = size)` (argument order changed) |
| `Cylinder(o, axis, r, ref)` crossing the cell | `Fibre(point, direction, r)` |
| `Ellipsoid(c, radii, [θ1, θ2], ref)` | `Rotate(Rotate(Ellipsoid(c, radii); axis = (0, 1, 0), angle = θ1); axis = (0, 0, 1), angle = θ2)` |
| `Fuse([a, b], ref)`, `Cut(a, b, ref)`, `Intersect(a, b, ref)` | `Fuse(a, b)`, `Cut(a, b)`, `Intersect(a, b)` |
| `createGmshModel` + `createMesh!` + `saveMesh(f)` | `generate(rve; output = f)` |
| `visualizeMesh()` | `generate(rve; view = true)` or `view_mesh(res)` |
| groups `Phase0`, `Phase<i>`, `Inclusions` | the matrix name and the phase names (no union group) |
| `Surfxmin`/`Linesxmin`/`Pointsxmin`, `Corners` | `xmin` (in dimensions 2, 1, 0), `corners` |
| (none) | `interface_<a>_<b>` (tags from 102) and `interfaces` (99): interfaces between phases |
| one group per inclusion | one group per phase (put inclusions in separate phases to keep them apart) |

Behaviour changes:

- Overlapping phases: the phase listed later wins (`overlaps = :paint`); `overlaps = :error`
  turns any overlap into an error.
- Mesh refinement around inclusions is opt-in per phase (`refine`); without it the mesh has
  the global size `h` (reduced only where curved boundaries need it). Refinement acts on the
  interfaces between phases, not on the cell faces, and grades with a constant ratio
  (`growth`); the size field is periodic.
- Periodic images are generated for every crossing of a periodic face, edge or corner,
  including inclusions larger than the cell.
- Tolerances are relative to the cell size: cells of any size and origin work.
- A mismatch between the meshes of opposite periodic faces is an error (`PeriodicityError`).

Example: the second example of the old README

```julia
# 2024
geo4 = Fuse([Cylinder([-0.1, 0.5, 0.5], [2.0, 0, 0], 0.1, 2),
             Cylinder([0.5, -0.1, 0.5], [0, 2.0, 0], 0.1, 2),
             Cylinder([0.5, 0.5, -0.1], [0, 0, 2.0], 0.1, 2)], 2)
geo7 = Intersect(Sphere([0.5, 0.5, 0.5], 0.25, 2), Box([L, L, L], [0.5 - L/2, 0.5 - L/2, 0.5 - L/2], 2), 2)
startGmsh()
model, tagvol = createGmshModel(RVE([1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1), (Cut(geo7, geo4, 2),), "RVE")
createMesh!(model, myRVE, myInclusions, tagvol); saveMesh("RVE.msh"); stopGmsh()

# version 0.2
c = (0.5, 0.5, 0.5)
cross = Fuse(Fibre(c, (1, 0, 0), 0.1), Fibre(c, (0, 1, 0), 0.1), Fibre(c, (0, 0, 1), 0.1))
core = Intersect(Sphere(c, 0.25), Box(center = c, size = (L, L, L)))
rve = RVE(Cell((1.0, 1.0, 1.0)), Phase(:inclusion, Cut(core, cross); refine = InterfaceRefinement(level = 2)))
generate(rve; h = 0.1, output = "RVE.msh")
```

Known limitations:

- Bodies flush with a periodic face: the part of their boundary lying on the face is treated
  as cell face (no interface refinement, no interface group there).
- With refined inclusions the periodic size field makes meshing 1.5–5 times slower
  (`MeshOptions(periodic_sizes = false)` to compare, a larger `growth` to speed it up).
- Volumes of pieces bounded by B-spline surfaces (ellipsoids, non-uniform scalings) cut by a
  face are computed by OpenCASCADE to about 1e-2.
