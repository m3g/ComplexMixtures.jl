import Pkg;
Pkg.activate(".");

using ComplexMixtures
using PDBTools

# Here we will project the contributions of the polymer atoms to the MDDF
# into the space around the polymer, and visualize the density of DMF

system = read_pdb("./equilibrated.pdb")
acr = select(system, "resname FACR or resname ACR or resname LACR")
results = load("./mddf.json")

# The PDB file does not contain the element column, and the element of
# the atoms is guessed from the atom names. Here, the CL (terminal methyl
# carbon) atom name would be interpreted as Chlorine. Thus we set the
# elements of the polymer atoms explicitly (all are C, H, N, or O):
for atom in acr
    atom.pdb_element = string(first(name(atom)))
end

# Compute the 3D density grid around the polymer
grid = grid3D(results, system; dmin=1.5, dmax=3.5)

# Regions associated with the hydrogen-bonding peak of the MDDF (d < 2.0 Å),
# with at least 10% of the maximum contribution. The polymer is shown in
# balls and sticks, and the grid points as dots colored by their relative
# contribution, from white (zero) to red (one)
hbonds = filter(p -> occup(p) < 2.0 && beta(p) > 0.1, grid)
view_hbonds = visualize(
    acr => (style=:ballandstick,),
    hbonds => (
        style=:dots, 
        color_by=beta.(hbonds), 
        colormap=:rwb, 
        color_range=(1, -1),
        opacity=0.8,
    ),
)
save("./grid_hbonds.html", view_hbonds)

# Regions associated with the second peak of the MDDF (2.0 Å < d < 3.5 Å),
# with at least 50% of the maximum contribution
second_peak = filter(p -> 2.0 < occup(p) < 3.5 && beta(p) > 0.5, grid)
view_second_peak = visualize(
    acr => (style=:ballandstick,),
    second_peak => (
        style=:dots, 
        color_by=beta.(second_peak), 
        colormap=:rwb, 
        color_range=(1, -1),
        opacity=0.8,
    ),
)
save("./grid_second_peak.html", view_second_peak)

# The grid can also be converted into volumetric data, smoothed with a
# Gaussian function of width 0.5 Å, and represented by isosurfaces. Here,
# the isosurfaces correspond to 50% (orange, transparent) and 75% (red) 
# of the maximum value of the smoothed data
density = volumetric_data(grid; sigma=0.5)
dmax = maximum(density.data)
view_isosurfaces = visualize(
    acr => (style=:ballandstick,),
    density => (isovalue=0.5 * dmax, color="orange", opacity=0.4),
    density => (isovalue=0.75 * dmax, color="red"),
)
save("./isosurfaces.html", view_isosurfaces)

# Write the volumetric data to a file, to be visualized in other software
write_dx("./density.dx", density)
println("Views saved to grid_hbonds.html, grid_second_peak.html, and isosurfaces.html")
