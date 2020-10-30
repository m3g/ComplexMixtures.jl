using Test
using Random
using PDBTools

using Revise
using ComplexMixtures
const CM = ComplexMixtures

dir="./test/data/NAMD"
atoms = readPDB("$dir/structure.pdb")  
options = CM.Options(lastframe=1,seed=1234567,nthreads=1,silent=true)

# Example 1: protein-tmao

# CM.save(R,"$dir/protein_tmao.json")
R_save = CM.load("$dir/protein_tmao.json")
protein = CM.Selection(select(atoms,"protein"),nmols=1)
tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
R = CM.mddf(traj,options)

#@test t = isapprox(R,R_save,debug=true) 

trajectory = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 

samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)
R = CM.Result(trajectory,options)
framedata = CM.FrameData(trajectory,R)
CM.nextframe!(trajectory)
iframe = 1
CM.mddf_frame!(iframe,framedata,options,R)





