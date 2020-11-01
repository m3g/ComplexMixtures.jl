using Revise
using ComplexMixtures, PDBTools
using Random
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

println(" --------------------------------------------------------------")
println(" Compiling - single thread ")
options = CM.Options(stride=1,seed=1234567,nthreads=1,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
@time R = CM.mddf(traj,options);
println(" Compiling - multi-threaded ")
options = CM.Options(stride=1,seed=1234567,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
@time R = CM.mddf(traj,options);
println(" Compiling - self single thread ")
options = CM.Options(stride=1,seed=1234567,nthreads=1,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",tmao) 
@time R = CM.mddf(traj,options);
println(" Compiling - self multi-threaded ")
options = CM.Options(stride=1,seed=1234567,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",tmao) 
@time R = CM.mddf(traj,options);

println(" --------------------------------------------------------------")
println(" Protein-TMAO - Single threaded ")
options = CM.Options(stride=1,seed=1234567,nthreads=1,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
@time R = CM.mddf(traj,options);


println(" --------------------------------------------------------------")
println(" Protein-TMAO - nthreads = $nthreads ")
options = CM.Options(stride=1,seed=1234567,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
@time R = CM.mddf(traj,options);
  
println(" --------------------------------------------------------------")
println(" Water-TMAO - Single threaded ")
options = CM.Options(stride=1,seed=1234567,nthreads=1,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",water,tmao) 
@time R = CM.mddf(traj,options);

println(" --------------------------------------------------------------")
println(" Water-TMAO - nthreads = $nthreads ")
options = CM.Options(stride=1,seed=1234567,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",water,tmao) 
@time R = CM.mddf(traj,options);
  
println(" --------------------------------------------------------------")
println(" Water-TMAO - Single threaded ")
options = CM.Options(stride=5,seed=1234567,nthreads=1,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",water) 
@time R = CM.mddf(traj,options);

println(" --------------------------------------------------------------")
println(" water-water - nthreads = $nthreads ")
options = CM.Options(stride=1,seed=1234567,silent=true)
traj = CM.Trajectory("$dir/trajectory.dcd",water) 
@time R = CM.mddf(traj,options);
println(" --------------------------------------------------------------")

  
