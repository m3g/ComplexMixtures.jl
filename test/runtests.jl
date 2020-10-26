
using Test
using ComplexMixtures, PDBTools
using Random
const CM = ComplexMixtures

# Example 1: protein-tmao with NAMD

dir="./data/NAMD"
R_save = CM.load("$dir/protein_tmao.json")
atoms = readPDB("$dir/structure.pdb","protein or resname TMAO")  
protein = CM.Selection(select(atoms,"protein"),nmols=1)
tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
options = CM.Options(stride=5)
R = CM.mddf(traj,options)
t = isapprox(R,R_save,debug=true) 
@test t


