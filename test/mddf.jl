@testset "NAMD" begin
  
   using PDBTools 
   atoms = PDBTools.readPDB("./data/NAMD/structure.pdb")
   solute = PDBTools.select(atoms,"protein")
   solvent = PDBTools.select(atoms,"resname TMAO")
   options = ComplexMixtures.Options(lastframe=1)
   trajectory = ComplexMixtures.Trajectory("./data/NAMD/trajectory.dcd",solute,solvent)
   R = ComplexMixtures.mddf(trajectory,options)

   @test  



end
