#
# Protein - Water
#

using PDBTools
include("../../src/ComplexMixtures.jl")

dir="/home/leandro/Documents/Emilia/corrected_concentration/50"
pdb="$dir/SC9_50%gly_corrected.pdb"
traj="$dir/01/production.dcd"

atoms = PDBTools.readPDB(pdb)
solute_indexes = PDBTools.select(atoms,"protein")
solute = ComplexMixtures.Selection( solute_indexes, nmols=1 )
solvent_indexes = PDBTools.select(atoms,"resname TIP3")
solvent = ComplexMixtures.Selection( solvent_indexes, natomspermol=3 )
trajectory = ComplexMixtures.Trajectory(traj,solute,solvent)
options = ComplexMixtures.Options(binstep=0.2,n_random_samples=10)
@time R = ComplexMixtures.mddf(trajectory,options)
