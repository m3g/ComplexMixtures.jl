
# The solute is a single protein molecule (infinte dilution case). In this case,
# use the option nmols=1
solute = MDDF.Solute( MDDF.VMDselect("structure.pdb","protein",vmd="/usr/local/bin/vmd"), 
                      nmols=1 )

# The solvent is TMAO, which has 14 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent = MDDF.Solvent( MDDF.VMDselect("structure.pdb","resname TMAO",vmd="/usr/local/bin/vmd"), 
                        natomspermol=14 ) 

# Initialize trajectroy data structure and open input stream
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)

# Run MDDF calculation 
#mddf = MDDF.mddf_naive(solute,solvent, 




