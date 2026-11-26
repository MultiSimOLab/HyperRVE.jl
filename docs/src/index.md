```@meta
CurrentModule = RVEs
```

# RVEs.jl

Periodic representative volume element (RVE) meshes with the Gmsh Julia API: describe a
[`Cell`](@ref), fill it with [`Phase`](@ref)s built from primitives, booleans and transforms,
and call [`generate`](@ref).

```julia
using RVEs
cell = Cell((1.0, 1.0, 1.0))
rve  = RVE(cell, Phase(:particles, [Sphere((0.5, 0.5, 0.5), 0.25), Sphere((0, 0, 0), 0.25)]; refine = 0.03))
res  = generate(rve; h = 0.1, output = "bcc.msh")
volume_fraction(res, :particles)          # ≈ 0.131
```

- cells of any size and origin, 2D and 3D, mixed periodicity, supercells;
- periodic images of inclusions crossing faces, edges and corners, items larger than the cell;
- phases painted in order, porous media without matrix;
- interface refinement with a periodic size field, phase sizes, order-2 elements, quality report;
- physical groups for phases, cell faces/edges/corners and interfaces, ready for GridapGmsh.

See the [User guide](guide.md), the [Examples](examples.md) and the
[migration notes](migration.md) for scripts written for the 2024 API.

## API

```@index
```

```@autodocs
Modules = [RVEs]
```
