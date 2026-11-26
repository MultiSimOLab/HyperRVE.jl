# Examples

The folder `examples/` of the repository contains one script per capability; each writes its
mesh next to the script.

| script | shows |
|--------|-------|
| `01_bcc_spheres.jl` | body-centred cubic spheres, periodic images, interface refinement |
| `02_fibres.jl` | hexagonal packing of infinite fibres, named matrix |
| `03_non_unit_cell.jl` | cell of any size and origin, mixed periodicity, rotated ellipsoid |
| `04_supercell.jl` | 3 × 3 × 3 supercell |
| `05_two_dimensional.jl` | 2D cell with disks and a rotated ellipse |
| `06_booleans.jl` | inclusion made with `Intersect`, `Cut` and `Fuse` |
| `07_mesh_control.jl` | phase sizes, graded refinement near a face, order 2, quality report |
| `08_coated_particle_interfaces.jl` | overlapping phases and interface groups |
| `09_porous.jl` | no matrix: porous solid and its pore walls |
| `10_strut_lattices.jl` | strut lattices: presets by density, in a matrix, custom topology, finite block |

Run them with `julia --project=. examples/01_bcc_spheres.jl` from the repository root.
