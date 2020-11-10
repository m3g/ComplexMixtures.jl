using Test
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

@testset "PDB" begin

  #
  # Tests with trajectory given in a PDB file
  #
  
  dir="./data/PDB"
  atoms = readPDB("$dir/trajectory.pdb","model 1")  
  options = CM.Options(stride=1,seed=321,StableRNG=true,nthreads=1,silent=true)
 
  # Example 1: protein-tmao
  
  # CM.save(R,"$dir/protein_tmao.json")
  R_save = CM.load("$dir/protein_tmao.json")
  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.pdb",protein,tmao,format="PDBTraj") 
  R = CM.mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 
  
  # Example 3: tmao-tmao
  
  # CM.save(R,"$dir/tmao_tmao.json")
  R_save = CM.load("$dir/tmao_tmao.json")
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.pdb",tmao,format="PDBTraj") 
  R = CM.mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 
  
end 


