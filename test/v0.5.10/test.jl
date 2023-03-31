import Pkg;
Pkg.activate(".");

using ComplexMixtures
using PDBTools

dir = "$(@__DIR__)/../data/NAMD"
const CM = ComplexMixtures

nthreads = 6
if nthreads == 1
    error(" These tests have to be run with multiple threads as well. ")
end

atoms = readPDB("$dir/structure.pdb")
protein = Selection(select(atoms, "protein"), nmols = 1)
tmao = Selection(select(atoms, "resname TMAO"), natomspermol = 14)
water = Selection(select(atoms, "water"), natomspermol = 3)

options_single = Options(stride = 1, seed = 321, nthreads = 1, silent = true, lcell = 2)
options_multi =
    Options(stride = 1, seed = 321, silent = true, nthreads = nthreads, lcell = 2)

println(" --------------------------------------------------------------")
println(" Compiling - single thread ")
traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
@time R = mddf(traj, options_single)
println(" Compiling - multi-threaded ")
traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
@time R = mddf(traj, options_multi)
println(" Compiling - self single thread ")
traj = Trajectory("$dir/trajectory.dcd", tmao)
@time R = mddf(traj, options_single)
println(" Compiling - self multi-threaded ")
traj = Trajectory("$dir/trajectory.dcd", tmao)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" Protein-TMAO - Single threaded ")
traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" Protein-TMAO - nthreads = $nthreads ")
traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" Protein-Water - Single threaded ")
traj = Trajectory("$dir/trajectory.dcd", protein, water)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" Protein-Water - nthreads = $nthreads ")
traj = Trajectory("$dir/trajectory.dcd", protein, water)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" Water-TMAO - Single threaded ")
traj = Trajectory("$dir/trajectory.dcd", water, tmao)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" Water-TMAO - nthreads = $nthreads ")
traj = Trajectory("$dir/trajectory.dcd", water, tmao)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" TMAO-TMAO - Single threaded ")
traj = Trajectory("$dir/trajectory.dcd", tmao)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" TMAO-TMAO - nthreads = $nthreads ")
traj = Trajectory("$dir/trajectory.dcd", tmao)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" Water-Water - Single threaded ")
traj = Trajectory("$dir/trajectory.dcd", water)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" Water-Water - nthreads = $nthreads ")
traj = Trajectory("$dir/trajectory.dcd", water)
@time R = mddf(traj, options_multi)
println(" --------------------------------------------------------------")