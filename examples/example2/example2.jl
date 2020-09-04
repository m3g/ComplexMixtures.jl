#
# Protein - Water
#

using ComplexMixtures ; const CM = ComplexMixtures
using PDBTools
using Plots
ENV["GKSwstype"] = "nul"

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
atoms = readPDB("../../test/data/NAMD/structure.pdb")

# The solute is a single protein molecule (infinte dilution case). In this case,
# use the option nmols=1
solute_indexes = select(atoms,"protein")
solute = CM.Selection( solute_indexes, nmols=1 )

# The solvent is Water, which has 3 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent_indexes = select(atoms,"resname TIP3")
solvent = CM.Selection( solvent_indexes, natomspermol=3 )

# Initialize trajectroy data structure and open input stream
trajectory = CM.Trajectory("../../test/data/NAMD/trajectory.dcd",solute,solvent)

# Input options for the calcualtion
options = CM.Options(binstep=0.2)

# Run CM calculation, and get the resutls in the R structure
@time R = CM.mddf(trajectory,options)

CM.save(R,"example1.json")
CM.write(R,"example1.dat",solute,solvent)

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
savefig("./example2.pdf")


