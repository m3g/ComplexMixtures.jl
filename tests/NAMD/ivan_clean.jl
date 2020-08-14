#
# Protein - TMAO (compare new and old implementations)
#

include("../../src/ComplexMixtures.jl")
using PDBTools

dir="/home/leandro/Drive/Alunos/Ivan/PCCP_revision"

atoms = PDBTools.readPDB("$dir/6Mnative.pdb")

protein = PDBTools.select(atoms,"protein")
solute = ComplexMixtures.Selection(protein,nmols=1 )

water = PDBTools.select(atoms,"water")
solvent = ComplexMixtures.Selection(water,natomspermol=3)

options = ComplexMixtures.Options(n_random_samples=10,lastframe=200)

trajectory = ComplexMixtures.Trajectory("$dir/6Mnative.dcd",solute,solvent)

@time lcP = ComplexMixtures.mddf(trajectory,options)


