#
# Water - TMAO
#

using Plots
nogtk()

include("../src/MDDF.jl")

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
using PDBTools
atoms = PDBTools.readPDB("./structure.pdb")

# The solute is water, with 3 atoms
solute_indexes = [ atom.index for atom in filter( atom -> (atom.resname == "TIP3"), atoms ) ]
solute = MDDF.Solute( solute_indexes, natomspermol=3 )

# The solvent is TMAO, which has 14 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TMAO", atoms ) ]
solvent = MDDF.Solvent( solvent_indexes, natomspermol=14 )

# Input options for the calcualtion
options = MDDF.Options(output="example.dat",binstep=0.2,n_random_samples=1000)

nlabel="lc"
rlabel="lcP"

# Run MDDF calculation, and get the resutls in the R structure
println("$nlabel:")
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
@time N = MDDF.mddf_linkedcells(trajectory,options)

# Run MDDF calculation, and get the resutls in the R structure
println("$rlabel:")
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
@time R = MDDF.mddf_linkedcells_parallel(trajectory,options)

#println("naive:")
#trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
#@time naive = MDDF.mddf_naive(trajectory,options)

plot(layout=(8,1))

sp=1
plot!(ylabel="MDDF",subplot=sp)
plot!(N.d,N.mddf,subplot=sp,label="$nlabel - mddf")
plot!(N.d,N.rdf,subplot=sp,label="$nlabel - rdf")
scatter!(R.d,R.mddf,subplot=sp,label="$rlabel - mddf")
scatter!(R.d,R.rdf,subplot=sp,label="$rlabel - rdf")
plot!(legend=:topright,subplot=sp)

sp=2
plot!(ylabel="KB",subplot=sp)
plot!(N.d,N.kb,subplot=sp,label="$nlabel - mddf")
plot!(N.d,N.kb_rdf,subplot=sp,label="$nlabel - rdf")
scatter!(R.d,R.kb,subplot=sp,label="$rlabel - mddf")
scatter!(R.d,R.kb_rdf,subplot=sp,label="$rlabel - rdf")
plot!(legend=:topright,subplot=sp)

sp=3
plot!(ylabel="Count",subplot=sp)
plot!(N.d,N.md_count,subplot=sp,label="$nlabel",linewidth=3,alpha=0.5)
scatter!(R.d,R.md_count,subplot=sp,label="$rlabel")
plot!(N.d,N.md_count_random,subplot=sp,label="$nlabel - rand",linewidth=3,alpha=0.5)
scatter!(R.d,R.md_count_random,subplot=sp,label="$rlabel -rand")

sp=4
plot!(ylabel="Shell vol", subplot=sp)
plot!(N.d,N.volume.shell,subplot=sp,label="$nlabel")
scatter!(R.d,R.volume.shell,subplot=sp,label="$rlabel")

sp=5
plot!(ylabel="Sum MD", subplot=sp)
plot!(N.d,N.sum_md_count,subplot=sp,label="$nlabel - md")
scatter!(R.d,R.sum_md_count,subplot=sp,label="$rlabel - md")

sp=6
plot!(ylabel="Sum RAND", subplot=sp)
plot!(N.d,N.sum_md_count_random,subplot=sp,label="$nlabel - rand")
plot!(N.d,N.sum_rdf_count,subplot=sp,label="$nlabel - rdf")
scatter!(R.d,R.sum_md_count_random,subplot=sp,label="$rlabel - rand")
scatter!(R.d,R.sum_rdf_count,subplot=sp,label="$rlabel - rdf")
plot!(legend=:topleft,subplot=sp)

sp=7
plot!(ylabel="Atom contributions", subplot=sp)
plot!(N.d,N.solute_atom,subplot=sp,label="$nlabel")
scatter!(R.d,R.solute_atom,subplot=sp,label="$rlabel")

sp=8
plot!(ylabel="Count RDF",subplot=sp)
plot!(N.d,N.rdf_count,subplot=sp,label="$nlabel",linewidth=3,alpha=0.5)
scatter!(R.d,R.rdf_count,subplot=sp,label="$rlabel")
plot!(N.d,N.rdf_count_random,subplot=sp,label="$nlabel - rand",linewidth=3,alpha=0.5)
scatter!(R.d,R.rdf_count_random,subplot=sp,label="$rlabel -rand")

plot!(size=(800,1500))
savefig("./plots.pdf")






