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
println("Views saved to grid_hbonds.html and grid_second_peak.html")
