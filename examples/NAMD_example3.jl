
include("../src/MDDF.jl")

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
using PDBTools
atoms = PDBTools.readPDB("./structure.pdb")

# Self correlation, thus the solute and solvent indexes are identical

solute_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TIP3", atoms ) ]
solute = MDDF.Solute( solute_indexes, natomspermol=3 )

solvent_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TIP3", atoms ) ]
solvent = MDDF.Solvent( solvent_indexes, natomspermol=3 )

# Initialize trajectroy data structure and open input stream
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)

# Input options for the calcualtion
options = MDDF.Options(output="example.dat",binstep=0.2,lastframe=-1)

# Run MDDF calculation, and get the resutls in the R structure
#N = MDDF.mddf_naive_self(trajectory,options)

R = MDDF.mddf_linkedcells_self(trajectory,options)

using Plots
nogtk()

plot(layout=(2,1))

x = [0,10]
y = [1, 1]

sp=1
plot!(ylabel="MDDF or RDF",subplot=sp)
scatter!(R.d,R.mddf,subplot=sp,label="mddf")
scatter!(R.d,R.rdf,subplot=sp,label="rdf")
plot!(N.d,R.mddf,subplot=sp,label="mddf - naive")
plot!(N.d,R.rdf,subplot=sp,label="rdf - naive")
plot!(x,y)
plot!(legend=:topright,subplot=sp)

sp=2
plot!(ylabel="KB",subplot=sp)
scatter!(R.d,R.kb,subplot=sp,label="mddf")
scatter!(R.d,R.kb_rdf,subplot=sp,label="rdf")
plot!(N.d,R.kb,subplot=sp,label="mddf - naive")
plot!(N.d,R.kb_rdf,subplot=sp,label="rdf - naive")
plot!(legend=:topright,subplot=sp)

plot!(size=(600,800))
savefig("./plots.pdf")
