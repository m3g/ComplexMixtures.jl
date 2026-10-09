import Pkg;
Pkg.activate(".");
using PDBTools
using ComplexMixtures

# PDB file of the system simulated
atoms = read_pdb("./system.pdb")
protein = select(atoms, "protein")

# Load results of a ComplexMixtures run
results = load("./glyc50_results.json")

# Compute the 3D density grid. Here we use dmax=3.5 such that the
# grid is not too large. Provide a file name (e.g. "./grid.pdb") as
# the third argument to also write the grid to a PDB file.
grid = grid3D(results, atoms; dmin=1.5, dmax=3.5)

# Regions associated with the hydrogen-bonding peak of the MDDF (d < 2.0 Å),
# with at least 10% of the maximum contribution. The protein is shown as a 
# white cartoon, the grid points as dots colored by their relative contribution
# from white (zero) to red (≥ 0.4), and the GLU and ASP residues as green balls and sticks
hbonds = filter(p -> occup(p) < 2.0 && beta(p) > 0.1, grid)
view_hbonds = visualize(
    protein => (color="white",),
    hbonds => (style=:dots, color_by=beta.(hbonds), colormap=:rwb, color_range=(0.4, -0.4)),
    protein => (selection="resname GLU ASP", style=:ballandstick, color="green"),
)
save("./grid_hbonds.html", view_hbonds)

# Regions associated with the second peak of the MDDF (2.0 Å < d < 3.5 Å),
# with the THR, ASN, and TRP residues as balls and sticks
second_peak = filter(p -> 2.0 < occup(p) < 3.5 && beta(p) > 0.1, grid)
view_second_peak = visualize(
    protein => (color="white",),
    second_peak => (style=:dots, color_by=beta.(second_peak), colormap=:rwb, color_range=(0.4, -0.4)),
    protein => (selection="resname THR ASN TRP", style=:ballandstick),
)
save("./grid_second_peak.html", view_second_peak)
