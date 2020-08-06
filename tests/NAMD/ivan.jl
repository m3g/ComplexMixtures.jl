#
# Protein - TMAO (compare new and old implementations)
#

using MDDF
using PDBTools
using Plots

dir="/home/leandro/Drive/Alunos/Ivan/PCCP_revision"

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
atoms = PDBTools.readPDB("$dir/6Mnative.pdb")

# Set solute and solvent

protein = PDBTools.select(atoms,"protein")
solute = MDDF.Selection(protein,nmols=1 )

water = PDBTools.select(atoms,"water")
solvent = MDDF.Selection(water,natomspermol=3)

options = MDDF.Options()

# Initialize trajectroy data structure and open input stream
#R = MDDF.mddf(trajectory,options)
#MDDF.save(R,"ivan.json")

trajectory = MDDF.Trajectory("$dir/6Mnative.dcd",solute,solvent)
lc = MDDF.mddf_linkedcells(trajectory,options)

trajectory = MDDF.Trajectory("$dir/6Mnative.dcd",solute,solvent)
lcP = MDDF.mddf_linkedcells_parallel(trajectory,options)

plot(layout=(5,1))

sp=1
plot!(ylabel="mddf",sp=sp)
plot!(lc.d,lc.mddf,linewidth=2,label="lc",subplot=sp)
scatter!(lcP.d,lcP.mddf,linewidth=2,label="lcP",subplot=sp)

sp=2
plot!(ylabel="md_count",sp=sp)
plot!(lc.d,lc.md_count,linewidth=2,label="lc",subplot=sp)
scatter!(lcP.d,lcP.md_count,linewidth=2,label="lcP",subplot=sp)

sp=3
plot!(ylabel="md_count_random",sp=sp)
plot!(lc.d,lc.md_count_random,linewidth=2,label="lc",subplot=sp)
scatter!(lcP.d,lcP.md_count_random,linewidth=2,label="lcP",subplot=sp)

sp=4
plot!(ylabel="rdf_count",sp=sp)
plot!(lc.d,lc.rdf_count,linewidth=2,label="lc",subplot=sp)
scatter!(lcP.d,lcP.rdf_count,linewidth=2,label="lcP",subplot=sp)

sp=5
plot!(ylabel="rdf_count_random",sp=sp)
plot!(lc.d,lc.rdf_count_random,linewidth=2,label="lc",subplot=sp)
scatter!(lcP.d,lcP.rdf_count_random,linewidth=2,label="lcP",subplot=sp)

plot!(size=(400,800))
savefig("./ivan.pdf")
savefig("./plots.pdf")
