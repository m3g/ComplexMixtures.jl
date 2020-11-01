using Test
using ComplexMixtures, PDBTools
using Random
const CM = ComplexMixtures

@testset "Merge" begin

  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  water = CM.Selection(select(atoms,"water"),natomspermol=3)
  
  # CM.save(R,"$dir/merged.json")
  R_save = CM.load("$dir/merged.json")

  options = CM.Options(firstframe=1,lastframe=2,seed=1234567,nthreads=1,silent=true)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao,water) 
  R1 = CM.mddf(traj,options)

  options = CM.Options(firstframe=3,lastframe=6,seed=1234567,nthreads=1,silent=true)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao,water) 
  R2 = CM.mddf(traj,options)

  R = CM.merge([R1,R2])
  t = isapprox(R,R_save,debug=true) 
  @test t

end

