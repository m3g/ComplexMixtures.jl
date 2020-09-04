#
# Self - Water
#

using ComplexMixtures ; const CM = ComplexMixtures
using PDBTools
using Plots
ENV["GKSwstype"] = "nul"

atoms = readPDB("../../test/data/NAMD/structure.pdb")

# Self correlation, thus the solute and solvent indexes are identical
water_atoms = select(atoms,"water")
water = CM.Selection( water_atoms, natomspermol=3)

# Input options for the calcualtion
options = CM.Options(binstep=0.2)

# Run MDDF calculation, and get the resutls in the R structure
nlabel="lcP"
trajectory = CM.Trajectory("../../test/data/NAMD/trajectory.dcd",water)
@time R = CM.mddf(trajectory,options)

plot(layout=(4,1))

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

sp=3         
plot!(ylabel="atom contrib",subplot=sp)
for i in 1:water.natomspermol
  plot!(R.d,R.solute_atom[:,i],subplot=sp,label="")
end
y1 = similar(R.d)
for i in 1:size(R.d,1)
  y1[i] = sum(R.solvent_atom[i,:])
end
plot!(R.d,y1,subplot=sp,label="sum")
plot!(legend=:topright,subplot=sp)

sp=4
plot!(ylabel="count",subplot=sp)
plot!(R.d,R.md_count,subplot=sp,label="R md")
plot!(legend=:topleft,subplot=sp)

plot!(size=(600,800))
savefig("./example3.pdf")
