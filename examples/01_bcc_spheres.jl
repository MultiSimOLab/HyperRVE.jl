# Body-centred cubic arrangement of spheres in a periodic unit cell, with interface refinement.
using RVEs

cell = Cell((1.0, 1.0, 1.0))
r = 0.25
spheres = [Sphere(cellpoint(cell, (0.5, 0.5, 0.5)), r),   # body centre
           Sphere(cellpoint(cell, (0.0, 0.0, 0.0)), r)]   # corner: its 7 images are added automatically
rve = RVE(cell, Phase(:particles, spheres; refine = 0.03))

res = generate(rve; h = 0.1, output = joinpath(@__DIR__, "bcc.msh"))
display(res)
println("\nvolume fraction of particles: ", volume_fraction(res, :particles),
        " (exact ", 2 * 4π * r^3 / 3, ")")
