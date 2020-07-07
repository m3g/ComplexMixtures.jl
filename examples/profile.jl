
using Profile

include("../src/MDDF.jl")

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
using PDBTools
atoms = PDBTools.readPDB("./structure.pdb")

# The solute is a single protein molecule (infinte dilution case). In this case,
# use the option nmols=1
solute_indexes = [ atom.index for atom in filter( atom -> PDBTools.isprotein(atom), atoms ) ]
solute = MDDF.Solute( solute_indexes, nmols=1 )

# Input options for the calcualtion
options = MDDF.Options(output="example.dat",binstep=0.2,lastframe=1)

Profile.clear()

for num in [ 10000000 ]

  println(" n = ", num) 

  solvent_indexes = [ atom.index for atom in filter( atom -> (atom.resname == "TIP3" && atom.resnum <= num ), atoms ) ]
  solvent = MDDF.Solvent( solvent_indexes, natomspermol=3 )

  # Initialize trajectroy data structure and open input stream

  # Run MDDF calculation, and get the resutls in the R structure
  println("linkedcells:")
  trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
  R = MDDF.mddf_linkedcells(trajectory,options)
  trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
  @time R = MDDF.mddf_linkedcells(trajectory,options)

  #println("naive:")
  trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
  R = MDDF.mddf_naive(trajectory,options)
  trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
  @time R = MDDF.mddf_naive(trajectory,options)

  println(" profiling linkedcells... ")
  trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
  #@profile R = MDDF.mddf_linkedcells(trajectory,options) 

  trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
  println(" profiling naive... ")
  @profile R = MDDF.mddf_naive(trajectory,options) 

end


