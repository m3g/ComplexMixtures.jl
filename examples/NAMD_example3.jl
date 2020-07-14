
include("../src/MDDF.jl")

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
using PDBTools
atoms = PDBTools.readPDB("./structure.pdb")

# Self correlation, thus the solute and solvent indexes are identical

resname = "TMAO"
natomspermol = 14
#resname = "TIP3"
#natomspermol = 3

solute_indexes = [ atom.index for atom in filter( atom -> atom.resname == resname && atom.resnum in [1,2], atoms ) ]
solute = MDDF.Solute( solute_indexes, natomspermol=natomspermol )

solvent_indexes = [ atom.index for atom in filter( atom -> atom.resname == resname && atom.resnum in [1,2], atoms ) ]
solvent = MDDF.Solvent( solvent_indexes, natomspermol=natomspermol )

# Input options for the calcualtion
options = MDDF.Options(output="example.dat",binstep=0.2,lastframe=1,n_random_samples=10000)

# Run MDDF calculation, and get the resutls in the R structure
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
N = MDDF.mddf_naive_self(trajectory,options)

trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
R = MDDF.mddf_linkedcells_self(trajectory,options)

using Plots
nogtk()

plot(layout=(4,1))

x = [0,10]
y = [1, 1]

sp=1
plot!(ylabel="MDDF or RDF",subplot=sp)
scatter!(R.d,R.mddf,subplot=sp,label="mddf")
scatter!(R.d,R.rdf,subplot=sp,label="rdf")
plot!(N.d,N.mddf,subplot=sp,label="mddf - naive")
plot!(N.d,N.rdf,subplot=sp,label="rdf - naive")
plot!(x,y)
plot!(legend=:topright,subplot=sp)

sp=2
plot!(ylabel="KB",subplot=sp)
scatter!(R.d,R.kb,subplot=sp,label="mddf")
scatter!(R.d,R.kb_rdf,subplot=sp,label="rdf")
plot!(N.d,N.kb,subplot=sp,label="mddf - naive")
plot!(N.d,N.kb_rdf,subplot=sp,label="rdf - naive")
plot!(legend=:topright,subplot=sp)

sp=3         
plot!(ylabel="atom contrib",subplot=sp)
i = 5
plot!(R.d,R.solute_atom[:,i],subplot=sp,label="new",linewidth=2)
plot!(R.d,R.solvent_atom[:,i],subplot=sp,label="new",linewidth=2)
scatter!(N.d,N.solute_atom[:,i],subplot=sp,label="naive")
scatter!(N.d,N.solvent_atom[:,i],subplot=sp,label="naive")
y1 = similar(R.d)
y2 = similar(R.d)
for i in 1:size(R.d,1)
  y1[i] = sum(R.solvent_atom[i,:])
  y2[i] = sum(N.solvent_atom[i,:])
end
plot!(R.d,y1,subplot=sp,label="sum")
scatter!(R.d,y2,subplot=sp,label="sum")
plot!(legend=:topright,subplot=sp)

sp=4
plot!(ylabel="count",subplot=sp)
scatter!(R.d,R.md_count,subplot=sp,label="R md")
scatter!(R.d,R.md_count_random,subplot=sp,label="R rand")
plot!(N.d,N.md_count,subplot=sp,label="N md")
plot!(N.d,N.md_count_random,subplot=sp,label="N rand")
plot!(legend=:topright,subplot=sp)

plot!(size=(600,800))
savefig("./plots.pdf")
