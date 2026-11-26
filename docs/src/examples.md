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
| `11_random_packings.jl` | random spheres (target fraction, size distribution), aligned fibres, disks |
| `12_voronoi_polycrystal.jl` | periodic Voronoi polycrystal, grain groups, 2D tessellation |
| `13_foams_and_voids.jl` | closed-cell foam at a density, open-cell foam, void phases |
| `14_particle_chains.jl` | field-structured MRE: particle chains, isotropic reference, variants, 2D |
| `15_more_lattices.jl` | BCCZ and FCCZ lattices; 2D square, triangular, honeycomb, kagome and re-entrant (auxetic) lattices |
| `16_arrays_and_laminates.jl` | regular particle arrays (SC, BCC, FCC, hexagonal 2D) and laminates |
| `17_auxetic_cells.jl` | chiral and anti-chiral lattices, rotating squares and triangles (auxetic 2D) |
| `18_particle_mixtures.jl` | hard + soft particle mixtures, chains with filler, clustered disks |

Run them with `julia --project=. examples/01_bcc_spheres.jl` from the repository root.
