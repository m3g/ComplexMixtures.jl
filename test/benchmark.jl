using Revise
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

nthreads = Threads.nthreads()
if nthreads == 1
  error(" These tests have to be run with multiple threads as well. ")
end

dir="$(@__DIR__ )/data/NAMD"
atoms = readPDB("$dir/structure.pdb")  
protein = CM.Selection(select(atoms,"protein"),nmols=1)
tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
water = CM.Selection(select(atoms,"water"),natomspermol=3)

options_single = CM.Options(stride=1,seed=321,nthreads=1,silent=true,lcell=2)
options_multi = CM.Options(stride=1,seed=321,silent=true,lcell=2)

println(" --------------------------------------------------------------")
println(" Compiling - single thread ")
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
@time R = CM.mddf(traj,options_single);
println(" Compiling - multi-threaded ")
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
@time R = CM.mddf(traj,options_multi);
println(" Compiling - self single thread ")
traj = CM.Trajectory("$dir/trajectory.dcd",tmao) 
@time R = CM.mddf(traj,options_single);
println(" Compiling - self multi-threaded ")
traj = CM.Trajectory("$dir/trajectory.dcd",tmao) 
@time R = CM.mddf(traj,options_multi);

println(" --------------------------------------------------------------")
println(" Protein-TMAO - Single threaded ")
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
@time R = CM.mddf(traj,options_single);

println(" --------------------------------------------------------------")
println(" Protein-TMAO - nthreads = $nthreads ")
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
@time R = CM.mddf(traj,options_multi);
  
println(" --------------------------------------------------------------")
println(" Protein-Water - Single threaded ")
traj = CM.Trajectory("$dir/trajectory.dcd",protein,water) 
@time R = CM.mddf(traj,options_single);

println(" --------------------------------------------------------------")
println(" Protein-Water - nthreads = $nthreads ")
traj = CM.Trajectory("$dir/trajectory.dcd",protein,water) 
@time R = CM.mddf(traj,options_multi);
  
println(" --------------------------------------------------------------")
println(" Water-TMAO - Single threaded ")
traj = CM.Trajectory("$dir/trajectory.dcd",water,tmao) 
@time R = CM.mddf(traj,options_single);

println(" --------------------------------------------------------------")
println(" Water-TMAO - nthreads = $nthreads ")
traj = CM.Trajectory("$dir/trajectory.dcd",water,tmao) 
@time R = CM.mddf(traj,options_multi);
  
println(" --------------------------------------------------------------")
println(" Water-Water - Single threaded ")
traj = CM.Trajectory("$dir/trajectory.dcd",water) 
@time R = CM.mddf(traj,options_single);

println(" --------------------------------------------------------------")
println(" Water-Water - nthreads = $nthreads ")
traj = CM.Trajectory("$dir/trajectory.dcd",water) 
@time R = CM.mddf(traj,options_multi);
println(" --------------------------------------------------------------")

  
