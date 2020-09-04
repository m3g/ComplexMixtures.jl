#
# Protein - TMAO (compare naive and linkedcells)
#

using ComplexMixtures ; const CM = ComplexMixtures
using PDBTools
using DelimitedFiles

using Plots
ENV["GKSwstype"] = "nul"

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
atoms = readPDB("../../test/data/NAMD/structure.pdb")

# The solute is a single protein molecule (infinte dilution case). In this case,
# use the option nmols=1
solute_indexes = select(atoms,"protein")
solute = CM.Selection( solute_indexes, nmols=1 )

# The solvent is TMAO, which has 14 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent_indexes = select(atoms,"resname TMAO")
solvent = CM.Selection( solvent_indexes, natomspermol=14 )

# Input options for the calcualtion
options = CM.Options(binstep=0.2)

# Run MDDF calculation, and get the resutls in the R structure
trajectory = CM.Trajectory("../../test/data/NAMD/trajectory.dcd",solute,solvent)
@time R = CM.mddf(trajectory,options)

plot(layout=(6,1))

sp=1
plot!(ylabel="MDDF",subplot=sp)
plot!(R.d,R.mddf,subplot=sp,label="new - mddf")
plot!(R.d,R.rdf,subplot=sp,label="new - rdf")
plot!(legend=:topright,subplot=sp)

sp=2
plot!(ylabel="KB",subplot=sp)
plot!(R.d,R.kb,subplot=sp,label="new - mddf")
plot!(R.d,R.kb_rdf,subplot=sp,label="new - rdf")
plot!(legend=:topright,subplot=sp)

sp=3
plot!(ylabel="Count",subplot=sp)
scatter!(R.d,R.md_count,subplot=sp,label="new")
scatter!(R.d,R.md_count_random,subplot=sp,label="new -rand")

sp=4
plot!(ylabel="Shell vol", subplot=sp)
plot!(R.d,R.volume.shell,subplot=sp,label="new")

sp=5
plot!(ylabel="Sum MD", subplot=sp)
scatter!(R.d,R.sum_md_count,subplot=sp,label="new - md")

sp=6
plot!(ylabel="Sum RAND", subplot=sp)
scatter!(R.d,R.sum_md_count_random,subplot=sp,label="new - rand")
scatter!(R.d,R.sum_rdf_count,subplot=sp,label="new - rdf")
plot!(legend=:topleft,subplot=sp)

plot!(size=(800,1300))
savefig("./example4.pdf")



