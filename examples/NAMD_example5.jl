#
# Protein - Water
#

using Plots
nogtk()

include("../src/MDDF.jl")

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
using PDBTools
atoms = PDBTools.readPDB("./structure.pdb")

# The solute is a single protein molecule (infinte dilution case). In this case,
# use the option nmols=1
solute_indexes = [ atom.index for atom in filter( atom -> PDBTools.isprotein(atom), atoms ) ]
solute = MDDF.Solute( solute_indexes, nmols=1 )

# The solvent is TMAO, which has 14 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TIP3", atoms ) ]
solvent = MDDF.Solvent( solvent_indexes, natomspermol=3 )

# Initialize trajectroy data structure and open input stream
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)

# Input options for the calcualtion
options = MDDF.Options(output="example.dat",binstep=0.2)

# Run MDDF calculation, and get the resutls in the R structure
@time R = MDDF.mddf_linkedcells(trajectory,options)

plot(layout=(4,1))

x = [0,10]
y = [1, 1]

sp=1
plot!(ylabel="MDDF or RDF",subplot=sp)
plot!(R.d,R.mddf,subplot=sp,label="mddf")
plot!(R.d,R.rdf,subplot=sp,label="rdf")
plot!(x,y)
plot!(legend=:topright,subplot=sp)

sp=2
plot!(ylabel="KB",subplot=sp)
plot!(R.d,R.kb,subplot=sp,label="mddf")
plot!(R.d,R.kb_rdf,subplot=sp,label="rdf")
plot!(legend=:topright,subplot=sp)

sp=3         
plot!(ylabel="atom contrib",subplot=sp)
for i in 1:solvent.natomspermol
  plot!(R.d,R.solvent_atom[:,i],subplot=sp,label="",linewidth=2)
end
plot!(legend=:topright,subplot=sp)

sp=4
plot!(ylabel="count",subplot=sp)
plot!(R.d,R.md_count,subplot=sp,label="R md")
plot!(R.d,R.md_count_random,subplot=sp,label="R rand")
plot!(legend=:topleft,subplot=sp)

plot!(size=(600,800))
savefig("./plots.pdf")


