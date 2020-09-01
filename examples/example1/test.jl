#
# Protein - TMAO (compare new and old implementations)
#

include("../../src/ComplexMixtures.jl")
using PDBTools
using Plots
using DelimitedFiles
ENV["GKSwstype"] = "nul" 

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
atoms = PDBTools.readPDB("../../test/data/NAMD/structure.pdb")

# The solute is a single protein molecule (infinte dilution case). In this case,
# use the option nmols=1
protein = PDBTools.select(atoms,"protein")
solute = ComplexMixtures.Selection(protein,nmols=1)

# The solvent is TMAO, which has 14 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
tmao = PDBTools.select(atoms,"resname TMAO")
solvent = ComplexMixtures.Selection(tmao,natomspermol=14)

# Initialize trajectroy data structure and open input stream
trajectory = ComplexMixtures.Trajectory("../../test/data/NAMD/trajectory.dcd",solute,solvent)

# Input options for the calculation
options = ComplexMixtures.Options(binstep=0.2)

# Run ComplexMixtures calculation, and get the results in the R structure
R = ComplexMixtures.mddf(trajectory,options)

# Save data for future loeading
ComplexMixtures.save(R,"example1.json")

# Save data in human-readable format
ComplexMixtures.write(R,"example1.dat",solute,solvent)

old = readdlm("./test_reference.dat",comments=true,comment_char='#')

plot(layout=(6,1))

sp=1
plot!(ylabel="mddf",subplot=sp)
plot!(old[:,1],old[:,2],subplot=sp,label="old")
plot!(R.d,R.mddf,subplot=sp,label="new - mddf")
plot!(R.d,R.rdf,subplot=sp,label="new - rdf")
plot!(legend=:topright,subplot=sp)

sp=2
plot!(ylabel="KB",subplot=sp)
plot!(old[:,1],old[:,3],subplot=sp,label="old")
plot!(R.d,R.kb,subplot=sp,label="new - mddf")
plot!(R.d,R.kb_rdf,subplot=sp,label="new - rdf")
plot!(legend=:topright,subplot=sp)

sp=3
plot!(ylabel="Count",subplot=sp)
plot!(old[:,1],old[:,4],subplot=sp,label="old")
scatter!(R.d,R.md_count,subplot=sp,label="new")
plot!(old[:,1],old[:,5],subplot=sp,label="old - rand")
scatter!(R.d,R.md_count_random,subplot=sp,label="new -rand")

sp=4
plot!(ylabel="Shell vol", subplot=sp)
plot!(old[:,1],old[:,8],subplot=sp,label="old")
plot!(R.d,R.volume.shell,subplot=sp,label="new")

sp=5
plot!(ylabel="Sum MD", subplot=sp)
plot!(old[:,1],old[:,6],subplot=sp,label="old - md")
scatter!(R.d,R.sum_md_count,subplot=sp,label="new - md")

sp=6
plot!(ylabel="Sum RAND", subplot=sp)
plot!(old[:,1],old[:,7],subplot=sp,label="old - rand")
scatter!(R.d,R.sum_md_count_random,subplot=sp,label="new - rand")
scatter!(R.d,R.sum_rdf_count,subplot=sp,label="new - rdf")
plot!(legend=:topleft,subplot=sp)

#sp=7
#plot!(ylabel="Count",subplot=sp)
#z = zeros(R.nbins)
#@. z = old[:,5]/R.md_count_random
#plot!(R.d,z,subplot=sp,label="old/rand")
#println(z[R.nbins])

plot!(size=(800,1300))
savefig("./test.pdf")


