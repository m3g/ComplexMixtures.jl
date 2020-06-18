
include("../src/MDDF.jl")

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
using PDBTools
atoms = PDBTools.readPDB("./structure.pdb")

# The solute is water, with 3 atoms
solute_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TIP3", atoms ) ]
solute = MDDF.Solute( solute_indexes, natomspermol=3 )

# The solvent is TMAO, which has 14 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TMAO", atoms ) ]
solvent = MDDF.Solvent( solvent_indexes, natomspermol=14 )

# Initialize trajectroy data structure and open input stream
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)

# Input options for the calcualtion
options = MDDF.Options(output="example.dat",n_random_samples=1)

# Run MDDF calculation, and get the resutls in the R structure
R = MDDF.mddf_naive(trajectory,options)


using Plots
nogtk()

plot(layout=(2,1))

sp=1
plot!(ylabel="MDDF or RDF",subplot=sp)
plot!(R.d,R.mddf,subplot=sp,label="mddf")
plot!(R.d,R.rdf,subplot=sp,label="rdf")
plot!(legend=:topright,subplot=sp)

sp=2
plot!(ylabel="KB",subplot=sp)
plot!(R.d,R.kb,subplot=sp,label="mddf")
plot!(R.d,R.kb_rdf,subplot=sp,label="rdf")
plot!(legend=:topright,subplot=sp)

plot!(size=(600,800))
savefig("./plots.pdf")
