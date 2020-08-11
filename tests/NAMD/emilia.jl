#
# Protein - Water
#

using PDBTools
include("../../src/MDDF.jl")

dir="/home/leandro/Documents/Emilia/corrected_concentration/50"
pdb="$dir/SC9_50%gly_corrected.pdb"
traj="$dir/01/production.dcd"

atoms = PDBTools.readPDB(pdb)
solute_indexes = PDBTools.select(atoms,"protein")
solute = MDDF.Selection( solute_indexes, nmols=1 )
solvent_indexes = PDBTools.select(atoms,"resname TIP3")
solvent = MDDF.Selection( solvent_indexes, natomspermol=3 )
trajectory = MDDF.Trajectory(traj,solute,solvent)
options = MDDF.Options(binstep=0.2,n_random_samples=10)
@time R = MDDF.mddf(trajectory,options)
