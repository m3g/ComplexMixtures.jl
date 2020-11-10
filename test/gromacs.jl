using Test
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

@testset "Gromacs" begin

  #
  # Tests with Gromacs-XSC trajectory
  #
  
  dir="./data/Gromacs"
  atoms = readPDB("$dir/system.pdb")  
  options = CM.Options(stride=5,seed=321,StableRNG=true,nthreads=1,silent=true)
  
  # Example 1: protein-EMIM
  
  # CM.save(R,"$dir/protein_EMI.json")
  R_save = CM.load("$dir/protein_EMI.json")
  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  traj = CM.Trajectory("$dir/trajectory.xtc",protein,emi) 
  R = CM.mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 
  
  # Example 1: EMIM-DCA
  
  # CM.save(R,"$dir/EMI_DCA.json")
  R_save = CM.load("$dir/EMI_DCA.json")
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  dca = CM.Selection(select(atoms,"resname NC"),natomspermol=5)
  traj = CM.Trajectory("$dir/trajectory.xtc",emi,dca) 
  R = CM.mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 
  
  # Example 1: EMIM-EMIM
  
  # CM.save(R,"$dir/EMI_EMI.json")
  R_save = CM.load("$dir/EMI_EMI.json")
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  traj = CM.Trajectory("$dir/trajectory.xtc",emi) 
  R = CM.mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 
  
end


