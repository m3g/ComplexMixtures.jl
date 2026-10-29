using ComplexMixtures, PDBTools, LazyArtifacts
const CM = ComplexMixtures

atoms = read_pdb(ComplexMixtures.test_file("protein_tmao.pdb"))
protein = AtomSelection(select(atoms, "protein"), nmols=1)
tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
water = AtomSelection(select(atoms, "water"), natomspermol=3)

options_single = Options(stride=1, seed=321, nthreads=1, silent=true)
options_multi = Options(stride=1, seed=321, silent=true)

nthreads = options_multi.nthreads

println(" --------------------------------------------------------------")
println(" Compiling - single thread ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), protein, tmao)
@time R = mddf(traj, options_single)
println(" Compiling - multi-threaded ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), protein, tmao)
@time R = mddf(traj, options_multi)
println(" Compiling - self single thread ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), tmao)
@time R = mddf(traj, options_single)
println(" Compiling - self multi-threaded ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), tmao)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" Protein-TMAO - Single threaded ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), protein, tmao)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" Protein-TMAO - nthreads = $nthreads ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), protein, tmao)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" Protein-Water - Single threaded ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), protein, water)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" Protein-Water - nthreads = $nthreads ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), protein, water)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" Water-TMAO - Single threaded ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), water, tmao)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" Water-TMAO - nthreads = $nthreads ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), water, tmao)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" TMAO-TMAO - Single threaded ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), tmao)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" TMAO-TMAO - nthreads = $nthreads ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), tmao)
@time R = mddf(traj, options_multi)

println(" --------------------------------------------------------------")
println(" Water-Water - Single threaded ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), water)
@time R = mddf(traj, options_single)

println(" --------------------------------------------------------------")
println(" Water-Water - nthreads = $nthreads ")
traj = ComplexMixtures.Trajectory(ComplexMixtures.test_file("protein_tmao.dcd"), water)
@time R = mddf(traj, options_multi)
println(" --------------------------------------------------------------")
