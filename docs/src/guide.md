# User guide

```@meta
CurrentModule = HyperRVE
```

## Workflow

1. Describe the **cell** with [`Cell`](@ref): size, origin and periodic directions.
2. Describe the **phases** with [`Phase`](@ref): a name and one or more geometries built from
   primitives, booleans and transforms ([`Geometry`](@ref)).
3. Wrap them in an [`RVE`](@ref) (background `matrix`, overlap rule, supercell tiling).
4. Call [`generate`](@ref), which returns an [`RVEMesh`](@ref) and writes the mesh files.

```julia
using HyperRVE
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

Strut copies whose part inside the cell is covered by a node sphere are not built (the body is
unchanged, the fuse 20–50 % faster). For lattices whose nodes are buried in the struts (octet),
`node_radius = 1.42 * radius` is about four times faster at a volume change below 1e-6; for
bcc, kelvin and diamond a larger node sphere changes the joint and the density.

## More lattices, particle arrays and laminates

- `lattice(cell, :bccz; density)` and `:fccz`: BCC and FCC with struts along the z edges
  (exact `relative_density` and `strut_radius` like the other presets).
- 2D strut lattices: `lattice(cell, kind; thickness)` with `kind` = `:square`, `:triangular`,
  `:honeycomb`, `:kagome` or `:reentrant` (auxetic; `angle` of the inclined struts below the
  horizontal, default 30°, and `ratio` = vertical / inclined strut length, default 2). The
  cell must have the proportions of `lattice_cell(kind; a)`, which returns it for strut
  length `a`: `c = lattice_cell(:honeycomb; a = 0.2); lattice(c, :honeycomb; thickness = 0.02)`.
  `struts` also works in 2D cells (rectangles and disks; `thickness = 2 radius`).
- `particle_array(cell, kind; fraction | radius)`: `:sc`, `:bcc`, `:fcc` spheres or `:square`,
  `:hex` disks (the hexagonal cell is `Cell((a, √3 a))`).
- `laminate(cell; axis, fractions, names, offset)`: one phase per layer (repeated names form one
  phase), for `matrix = nothing`; by default the first layer is centred on the cell face so that
  no interface lies on a face.
- `relative_density(kind, t / a)` gives the exact area fraction of the 2D square, triangular,
  honeycomb and kagome lattices.

## Auxetic 2D cells

- `chiral(cell, kind; ring_radius, thickness, mirror = false, filled = false)`: rings of mean
  radius `ring_radius` joined by tangent ligaments, all of width `thickness`. `kind` is
  `:hexachiral`, `:tetrachiral`, `:trichiral` or the anti-chiral `:antitetrachiral`,
  `:antitrichiral`. `mirror` flips the handedness of the chiral kinds (the anti-chiral ones
  are their own mirror image) and `filled` makes the rings solid disks.
  The cell is `lattice_cell(kind; a)` with `a` the distance between connected rings.
- `rotating_squares(cell; angle, hinge)` and `rotating_triangles(cell; angle, hinge)`: rigid
  units of side `a` turned by `angle` (radians) from the closed state and joined at the
  corners by disks of diameter `hinge`; the cell is `lattice_cell(kind; a, angle)`. The cell
  proportions do not depend on `angle`, so `a` is deduced from the cell and `angle`: give the
  same `angle` to `lattice_cell` and to the function.
- `relative_density(kind; ring_radius, thickness, a)` and `relative_density(kind; angle, hinge, a)`
  give the exact solid fraction while neighbouring ligaments do not overlap
  (`thickness < 2 ring_radius / 3` for `:hexachiral`; with `filled`, `thickness <= 2 ring_radius`
  for the tetra- and `<= 6 ring_radius` for the tri- kinds); outside that range they throw.
- Nearly closed openings and very small hinges give sharp corners: refine with a smaller `h`
  or `MeshOptions(hmin, curvature)` if the element quality matters.

```julia
c = lattice_cell(:hexachiral; a = 1.0)
Phase(:lattice, chiral(c, :hexachiral; ring_radius = 0.2, thickness = 0.04); fuse = true)
c = lattice_cell(:rotating_squares; a = 1.0, angle = π / 6)
Phase(:squares, rotating_squares(c; angle = π / 6, hinge = 0.08); fuse = true)
```

Hollow struts: `hollow_lattice(cell, kind; radius, wall | density, node = :hollow)` (or
`hollow_struts(cell, nodes, edges; radius, wall)`) returns `(shell, core)`: the lattice of
outer `radius` and its cavity, to be painted over it by a later phase, void (empty tubes) or
material (core-shell). The cavity phase must come after the shell in the list of phases and
the RVE must keep the default `overlaps = :paint` (listed before, the lattice comes out solid).
`node = :solid` keeps the nodes solid, one closed cavity per strut (no `fuse` needed for it).
`relative_density(kind, x; wall, node)` is the exact density.

```julia
hl = hollow_lattice(cell, :octet; radius = 0.06, wall = 0.02)
RVE(cell, [Phase(:lattice, hl.shell; fuse = true), Phase(:cavity, hl.core; fuse = true, void = true)]; matrix = nothing)
```

## Plate lattices

`plate_cells(cell, kind; thickness | density, ratio)` returns the closed cells left by plates
on the planes of `kind`: `:sc` (cube faces), `:bcc` ({110}), `:fcc` ({111}) or a tuple such as
`(:sc, :fcc)`. The plates are the rest of the cell, so the cells go to a void phase and the
plates are the matrix (or the cells are a second material):

```julia
cells = plate_cells(cell, (:sc, :fcc); density = 0.3)          # isotropic SC + FCC by default
RVE(cell, Phase(:cells, cells; void = true); matrix = :plates)
```

`plate_density(cell, kind; thickness)` is the exact relative density and
`plate_thickness(cell, kind, density; ratio)` its inverse (`plate_density(kind, t / a)` and
`plate_thickness(kind, density)` for a cube of edge `a`). `ratio = :isotropic` gives the
elastically isotropic SC + FCC lattice (`t_sc / t_fcc = 8√3/9`); other combinations need
`ratio`. Thin plates need a mesh size of about half their thickness.

## Triply periodic minimal surfaces (experimental)

`tpms_level(kind; type, density, periods, resolution)` and `tpms_density(kind; type, level, ...)`
relate the level of a TPMS solid to its relative density: `kind` is `:gyroid`, `:schwarz_p`,
`:diamond`, `:iwp` or `:neovius`; `type = :network` is the solid `f ≤ level` and `:sheet` the
wall `|f| ≤ level`. The solid is the faceted surface of a periodic marching-tetrahedra grid with
`resolution` cubes per period (24 by default), smoothed onto the surface; the density is that of
this faceted solid. Networks below their pinch-off (Schwarz P below about 0.21) and walls thinner
than about two grid cubes raise an error.

```julia
tpms_level(:gyroid; type = :sheet, density = 0.3)      # ≈ 0.463
```

`tpms(kind; type, density | level, periods, resolution)` makes a TPMS phase item. It is meshed by
its own path (the triangulated surface is kept as it is and Gmsh fills the volumes), so it must
be the only item of the only phase, with a matrix or `matrix = nothing`, in a cubic periodic cell,
without tiling, `h`, `refine`, `fuse`, `void` or order 2. The mesh has the phase groups, the
interface, the faces, edges and corners, and every node of opposite faces is paired in the
`$Periodic` section of the `.msh` file (GridapGmsh reads it).

```julia
rve = RVE(cell, Phase(:gyroid, tpms(:gyroid; type = :sheet, density = 0.3)))
res = generate(rve; output = "gyroid.msh")
```

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

Mixtures and clusters:

- `avoid = items` packs a new population around inclusions already placed (another packing,
  `particle_chains`, or geometries made by hand), at least `min_gap` from them, periodic images
  included. One call and one phase per population, the largest first; each `fraction` refers to
  the whole cell. Spheres, disks, capsules and axis-parallel fibres are avoided exactly;
  cylinders, ellipsoids and ellipses through the capsule that contains them.
- `cluster = (n = 4, radius = 0.2)` places the inclusions uniformly in balls of radius `radius`
  around `n` random centres (or `centers = [...]`); `share = s` sends only a fraction `s` of the
  attempts to the clusters.

```julia
hard = random_spheres(cell; radius = 0.14, fraction = 0.10, rng = 1)
soft = random_spheres(cell; radius = 0.09, fraction = 0.06, avoid = hard, rng = 2)
RVE(cell, [Phase(:hard, hard), Phase(:soft, soft)]; overlaps = :error)
random_disks(Cell((2.0, 2.0)); radius = 0.04, fraction = 0.15, cluster = (n = 5, radius = 0.3), rng = 5)
```

Short fibres: `random_capsules(cell; length, radius, fraction | n, orientation)` places
non-overlapping capsules (`length` between the cap centres) with `orientation = :isotropic`,
`:planar` (with `axis`), a direction, a 3×3 orientation tensor `A` (symmetric, positive
semidefinite, trace 1; directions from the central angular Gaussian with second moment `A`),
or a function `rng -> direction`; the other keywords as
above, `avoid` included (hybrid particle + fibre composites). `orientation_tensor(fibres)` gives
the second-order orientation tensor. Random sequential addition of slender isotropic fibres
saturates early (about 0.15 for a total aspect ratio of 10); aligned fibres reach about 0.3.
Refine the fibre phase to `radius / 2`–`radius`.

```julia
f = random_capsules(cell; length = 0.3, radius = 0.04, fraction = 0.025, rng = 1)
RVE(cell, Phase(:fibres, f; refine = 0.03))
orientation_tensor(f)          # ≈ I/3
```

## Particle chains (magnetorheological elastomers)

`particle_chains(cell; axis, radius, gap, fraction | n_chains)` places spheres (disks in 2D)
in chains along a cell axis, the microstructure of field-structured MREs:

```julia
cell   = Cell((1.0, 1.0, 1.0))
chains = particle_chains(cell; axis = 3, radius = 0.05, gap = 0.025, fraction = 0.15, rng = 1)
rve    = RVE(cell, Phase(:particles, chains; refine = 0.015))
```

- Along a periodic axis the chains are infinite: the gap is adjusted upwards so that the chain
  closes (warning when it grows by more than half).
- `arrangement = :random` (random positions of the chains in the cross-section), `:square` or
  `:hex` (with `n_chains = (n1, n2)`); `stagger = :alternate` shifts neighbouring chains by half
  a pitch (body-centred tetragonal).
- `waviness`, `occupancy < 1` (broken chains), `aspect` (spheroids elongated along the chain).
- Mesh: the ligament between two particles of a chain is `gap` wide. `MeshOptions(gaps = 2)`
  puts two elements across every thin gap between spheres or disks (and between them and the
  faces), with a size that follows the local width of the ligament, and keeps the rest coarse
  (about 400 r / gap tetrahedra per gap; `optimize = :netgen` removes the occasional sliver).
  Refining the whole phase (`refine = gap / 2`) also resolves the ligaments but costs 10–50
  times more elements. `GapRefinement(elements, faces, growth)` sets the details; other
  shapes (ellipsoids, capsules) are not refined yet.
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
