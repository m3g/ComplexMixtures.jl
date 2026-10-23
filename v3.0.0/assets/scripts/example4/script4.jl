import Pkg;
Pkg.activate(".");

using PDBTools
using ComplexMixtures

# Here we will project the contributions of the glycerol atoms to the
# MDDF of water relative to glycerol into the space around one
# glycerol molecule, and visualize the density of water

# Load the PDB file of the system and the previously computed results
system = read_pdb("./equilibrated.pdb")
mddf_glyc_water = load("./mddf_glyc_water.json")

# The solute (glycerol) is composed of many molecules: the grid is built
# around one of them (the first one, by default; it can be chosen with
# the `molecule` keyword argument of grid3D). The contributions of the
# atoms are averaged over all glycerol molecules.
grid = grid3D(mddf_glyc_water, system; dmin=1.5, dmax=3.5)
glycerol = system[mddf_glyc_water.solute.indices[1:14]] # first molecule

# Regions associated with the hydrogen-bonding peak of the MDDF (d < 2.0 Å),
# with at least 10% of the maximum contribution. The glycerol molecule is
# shown in balls and sticks, and the grid points as dots colored by their
# relative contribution, from white (zero) to red (one)
hbonds = filter(p -> occup(p) < 2.0 && beta(p) > 0.1, grid)
view_hbonds = visualize(
    glycerol => (style=:ballandstick,),
    hbonds => (
        style=:dots,
        color_by=beta.(hbonds),
        colormap=:rwb,
        color_range=(1, -1),
        opacity=0.8,
    ),
)
save("./grid_hbonds.html", view_hbonds)

# The grid is converted into volumetric data, smoothed with a Gaussian
# function of width 0.5 Å, and represented by isosurfaces. Here, the
# isosurfaces correspond to 50% (orange, transparent) and 75% (red) of
# the maximum value of the smoothed data
density = volumetric_data(grid; sigma=0.5)
dmax = maximum(density.data)
view_isosurfaces = visualize(
    glycerol => (style=:ballandstick,),
    density => (isovalue=0.5 * dmax, color="orange", opacity=0.4),
    density => (isovalue=0.75 * dmax, color="red"),
)
save("./isosurfaces.html", view_isosurfaces)

# Write the volumetric data to a file, to be visualized in other software
write_dx("./density.dx", density)
println("Views saved to grid_hbonds.html and isosurfaces.html")
