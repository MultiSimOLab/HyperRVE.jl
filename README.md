# RVEs.jl


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://multisimolab.github.io/RVEs.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://multisimolab.github.io/RVEs.jl/dev/)
[![Build Status](https://github.com/MultiSimOLab/RVEs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MultiSimOLab/RVEs.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Modelling of periodic RVEs using GMSH

**Gmsh:**
Gmsh stands as an open-source 3D finite element mesh generator equipped with an integrated CAD engine and post-processor. It's designed for speed, efficiency, and user-friendliness, featuring parametric input and adaptable visualization capabilities. The software architecture revolves around four core modules: geometry, mesh, solver, and post-processing. These modules offer flexible control options, including a graphical interface, command-line interface, dedicated .geo files, and support for APIs in C++, C, Python, Julia, and Fortran.

**RVEs Package:**
The RVEs package serves as a versatile tool for mesh generation, leveraging the Gmsh-Julia-API for seamless integration. It provides a robust framework for automating mesh generation for custom model types, enabling users to automate the creation of intricate structures such as representative volume elements (RVEs). RVEs simplifies the mesh modeling process into several key stages:

1. **Geometry Setup:** Define geometric shapes using fundamental entities and boolean operations.
2. **Integration with Gmsh:** Incorporate shapes into Gmsh, execute boolean operations, and define physical groups.
3. **Mesh Generation:** Generate meshes with customizable refinement fields.
4. **Export and Visualization:** Export meshes into various formats and visualize the results.

**Current Focus and Future Expansion:**
At present, RVEs excels in automating the generation of representative volume elements containing multiple inclusion objects. However, its capabilities can be extended by adding support for additional geometric shapes and model types. It's essential to clarify that RVEs does not aim to replace the Gmsh scripting language. Instead, it acts as a facilitator, simplifying the automation of complex meshing tasks within a user-friendly programming environment such as Julia.

# MiniApps included  

 <p align="center"> 
&nbsp; &nbsp; &nbsp; &nbsp;
<img alt="Dark"
src="https://github.com/MultiSimOLab/RVEs.jl/blob/main/docs/imgs/rves.png" width="100%">
</p>






# Quick start

```julia
using RVEs

cell = Cell((1.0, 1.0, 1.0))                     # periodic unit cube (any size and origin work)

# body-centred cubic arrangement: the corner sphere is completed by its periodic images
spheres = [Sphere((0.5, 0.5, 0.5), 0.25), Sphere((0.0, 0.0, 0.0), 0.25)]
rve = RVE(cell, Phase(:particles, spheres; refine = InterfaceRefinement(level = 2)))

res = generate(rve; h = 0.1, output = "bcc.msh")
volume_fraction(res, :particles)                 # 2 · 4π/3 · 0.25³ ≈ 0.131
```

Inclusions built with boolean operations and fibres crossing the whole cell:

```julia
using RVEs

cell  = Cell((1.0, 1.0, 1.0))
c     = cellpoint(cell, (0.5, 0.5, 0.5))         # fractional -> physical coordinates
L     = 0.38
cross = Fuse(Fibre(c, (1, 0, 0), 0.1), Fibre(c, (0, 1, 0), 0.1), Fibre(c, (0, 0, 1), 0.1))
core  = Intersect(Sphere(c, 0.25), Box(center = c, size = (L, L, L)))

rve = RVE(cell, [Phase(:fibre, Fibre(cellpoint(cell, (0, 0.5, 0)), (1, 0, 0), 0.1)),
                 Phase(:inclusion, Cut(core, cross); refine = InterfaceRefinement(level = 2))])
res = generate(rve; h = 0.1, output = "RVE.msh")
```

`generate` opens and closes Gmsh by itself and returns plain Julia data (phase volumes and
fractions, physical groups, mesh statistics, periodicity checks). Every mesh written is
periodic node by node on opposite faces and carries the physical groups

| name | dimensions | content |
|------|-----------|---------|
| phase names (`matrix`, ...) | 3 | volumes of each phase |
| `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax` | 2, 1, 0 | closed cell faces (tags 1–6) |
| `corners`, `boundary` | 0 / 2, 1, 0 | the 8 corners (7), the whole cell surface (8) |
| `edge_ymin_zmin`, ... | 1, 0 | the 12 cell edges |
| `corner_xmin_ymin_zmin`, ... | 0 | each corner |

Version 2 is a complete rewrite; see `docs/src/migration.md` to port scripts written for the
2024 API.

# Project funded by:

- Grant PID2022-141957OA-C22 funded by MCIU/AEI/ 10.13039/501100011033  and by ''ERDF A way of making Europe''

 <p align="center"> 
&nbsp; &nbsp; &nbsp; &nbsp;
<img alt="Dark"
src="https://github.com/MultiSimOLab/RVEs.jl/blob/main/docs/imgs/aei.png" width="60%">
</p>