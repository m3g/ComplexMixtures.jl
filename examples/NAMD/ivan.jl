#
# Protein - TMAO (compare new and old implementations)
#

using MDDF
using PDBTools
using Plots

dir="/home/leandro/Drive/Alunos/Ivan/PCCP_revision"

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
atoms = PDBTools.readPDB("$dir/6Mnative.pdb")

# Set solute and solvent

protein = PDBTools.select(atoms,"protein")
solute = MDDF.Selection(protein,nmols=1 )

water = PDBTools.select(atoms,"water")
solvent = MDDF.Selection(water,natomspermol=3)

# Initialize trajectroy data structure and open input stream
trajectory = MDDF.Trajectory("$dir/6Mnative.dcd",solute,solvent)

R = MDDF.mddf(trajectory)
MDDF.save(R,"ivan.json")

plot(R.d,R.mddf,linewidth=2,label="")
savefig("./plots.pdf")
