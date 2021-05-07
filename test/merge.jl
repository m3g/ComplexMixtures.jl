using Test
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

@testset "Merge" begin

  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  
  tmao = Selection(select(atoms,"resname TMAO"),natomspermol=14)
  water = Selection(select(atoms,"water"),natomspermol=3)
  
  # save(R,"$dir/merged.json")
  R_save = load("$dir/merged.json")

  options = Options(firstframe=1,lastframe=2,seed=321,StableRNG=true,nthreads=1,silent=true)
  traj = Trajectory("$dir/trajectory.dcd",tmao,water) 
  R1 = mddf(traj,options)

  options = Options(firstframe=3,lastframe=6,seed=321,StableRNG=true,nthreads=1,silent=true)
  traj = Trajectory("$dir/trajectory.dcd",tmao,water) 
  R2 = mddf(traj,options)

  R = merge([R1,R2])
  @test isapprox(R,R_save,debug=true) 

end

