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
| `Polyhedron(vertices, faces)` (planar faces, closed surface; prisms, grains, faceted particles) | `Polygon(vertices)` |
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

## Strut lattices

`lattice(cell, kind; radius | density)` returns the items of a preset strut lattice (`:cubic`,
`:bcc`, `:fcc`, `:octet`, `:kelvin`, `:diamond`): one `Cylinder` per strut and one `Sphere`
per node (`node_radius = :auto` uses the largest strut at the node; `nothing` removes them).
Put them in a fused phase:

```julia
cell = Cell((1.0, 1.0, 1.0))
rve  = RVE(cell, Phase(:lattice, lattice(cell, :octet; density = 0.15); fuse = true); matrix = nothing)
```

`density` is the relative density of the fused lattice (exact formula for the presets in a cubic
cell, with the node spheres; `strut_radius(kind, density)` and `relative_density(kind, x)`
convert, with `x = radius / edge`). Any topology works with `struts(cell, nodes, edges; radius)`:
nodes in fractional coordinates and struts `(i, j)` or `(i, j, offset)`, where `offset` is the
lattice shift of node `j`; `lattice_topology(kind)` returns those of a preset.
`relative_density(cell, nodes, edges; radius)` is a first-order estimate (it counts overlaps at
the nodes twice). In a non-periodic cell the lattice is the complete finite truss, outer faces
included. Struts of different radii: `radius = [r1, r2, ...]`, one per strut.

## Random packings

`random_spheres(cell; ...)`, `random_disks(cell; ...)` (2D) and `random_fibres(cell; axis, ...)`
return non-overlapping inclusions placed by random sequential addition, ready for a phase
(without `fuse`):

```julia
cell = Cell((1.0, 1.0, 1.0))
rve  = RVE(cell, Phase(:particles, random_spheres(cell; radius = 0.08, fraction = 0.2, rng = 1)))
```

Give two of `radius`, `fraction` and `n`, a vector of radii, or a sampler `rng -> radius` with
`fraction` or `n`. `min_gap` (default `0.1 *` the smallest radius) keeps the surfaces apart;
`face_gap` (default `min_gap`) rejects inclusions nearly tangent to a face, edge or corner of
the cell, or centred nearly on a face plane: these give slivers and Gmsh periodicity failures.
`boundary = :inside` keeps every inclusion inside the cell (leaving an inclusion-free layer
next to the faces); along non-periodic directions `:periodic` lets inclusions cross the face, where
they are cut. `volume_fraction(items, cell)` gives the fraction of the returned inclusions. `rng` is an `AbstractRNG` or an
integer seed. Random sequential addition saturates at about 0.38 for equal spheres and 0.55 for
equal disks (less with gaps: about 0.30 and 0.48 with the default `min_gap`); above that the
call fails (an `ErrorException`) with what it placed (`strict = false`: warning and partial
packing).
Dense packings leave ligaments of width `min_gap` between inclusions: the mesh must resolve
them (`refine` on the phase or a smaller `h`), or raise `min_gap`; otherwise Gmsh may leave
pieces unmeshed and `generate` reports it.

## Particle chains (magnetorheological elastomers)

`particle_chains(cell; axis, radius, gap, fraction | n_chains)` places spheres (disks in 2D)
in chains along a cell axis, the microstructure of field-structured MREs:

```julia
cell   = Cell((1.0, 1.0, 1.0))
chains = particle_chains(cell; axis = 3, radius = 0.06, gap = 0.012, fraction = 0.15, rng = 1)
rve    = RVE(cell, Phase(:particles, chains; refine = 0.012))
```

- Along a periodic axis the chains are infinite: the gap is adjusted upwards so that the chain
  closes (warning when it grows by more than half).
- `arrangement = :random` (random positions of the chains in the cross-section), `:square` or
  `:hex` (with `n_chains = (n1, n2)`); `stagger = :alternate` shifts neighbouring chains by half
  a pitch (body-centred tetragonal).
- `waviness`, `occupancy < 1` (broken chains), `aspect` (spheroids elongated along the chain).
- Mesh: the ligament between two particles of a chain is `gap` wide; refine the phase to about
  `gap / 2`…`gap` and keep `gap >= 0.1 radius`. Dense chains with thin gaps give large meshes and
  Gmsh may fail on them ("overlapping facets"): `gap = 0.25 radius` with `refine = 0.6 gap` is a
  robust starting point; a larger `chain_gap` keeps neighbouring chains apart.
- The isotropic reference at the same fraction is `random_spheres(cell; radius, n = length(chains))`.

## Voronoi polycrystals

`voronoi(cell; n, rng)` (or `seeds = [...]`) computes a periodic Voronoi tessellation in pure
Julia; `min_distance` keeps random seeds apart (default half the mean grain size, `0` for
Poisson–Voronoi) and `lloyd = k` makes the grain sizes more uniform. `grain_phases(t)` gives one
phase per grain (`grain_1`, …) made of a `Polyhedron` (2D: `Polygon`); use them without matrix:

```julia
t   = voronoi(Cell((1.0, 1.0, 1.0)); n = 20, rng = 1)
rve = RVE(t.cell, grain_phases(t); matrix = nothing)
```

The grain boundaries are the interface groups `interface_grain_i_grain_j`. `groups` assigns
grains to phases (`grain_phases(t; groups = [1, 2, 1, ...])`, fused per phase).
`grain_volumes(t)` and `grains(t)` give the exact volumes and the polyhedra. `refine` and `h`
are passed to every grain phase. `voronoi(cell; ..., periodic = true)` on a cell with free faces
gives a block cut from a periodic tessellation. Random tessellations have very short edges;
they mesh, but leave a few poor elements. Seeds on a regular lattice perturbed by less than
about 1e-2 of the spacing give grain faces nearly coplanar with the cell faces, which Gmsh may
fail to make periodic: use the exact lattice or a larger perturbation. Many grains give many
interface groups (about 7 per grain).

## Void phases and foams

`Phase(name, items; void = true)` makes empty space: what it covers is removed from the mesh
(also over earlier phases), and its surfaces are the groups `interface_<a>_void`.

From a tessellation `t = voronoi(cell; n, rng)`:
- closed cells: `pores(t; wall)` are the grains shrunk by `wall / 2`, and
  `wall_thickness(t, density)` the wall for a relative density (exact):
  `RVE(cell, Phase(:pores, pores(t; wall = w); void = true); matrix = :walls)`;
- open cells: `foam_struts(t; radius)` are the edges as struts (short edges merged), for a
  fused phase: `RVE(cell, Phase(:foam, foam_struts(t; radius = 0.02); fuse = true); matrix = nothing)`
  (fully periodic tessellations; practical for about 5–20 cells). The geometry is exact, but
  meshing open foams from random tessellations is fragile (very small pieces of struts cut by
  the cell faces): regular seeds (BCC → Kelvin foam) mesh well.

`generate(...; interfaces = :union)` creates only the union group `interfaces`, useful with
many grains.

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
