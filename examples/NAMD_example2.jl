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
options = MDDF.Options(output="example.dat",binstep=0.2)

# Run MDDF calculation, and get the resutls in the R structure
println("Naive:")
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
@time N = MDDF.mddf_naive(trajectory,options)

# Run MDDF calculation, and get the resutls in the R structure
println("LinkedCells:")
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)
@time R = MDDF.mddf_linkedcells(trajectory,options)

plot(layout=(7,1))

sp=1
plot!(ylabel="MDDF",subplot=sp)
plot!(N.d,N.mddf,subplot=sp,label="naive - mddf")
plot!(N.d,N.rdf,subplot=sp,label="naive - rdf")
plot!(R.d,R.mddf,subplot=sp,label="new - mddf")
plot!(R.d,R.rdf,subplot=sp,label="new - rdf")
plot!(legend=:topright,subplot=sp)

sp=2
plot!(ylabel="KB",subplot=sp)
plot!(N.d,N.kb,subplot=sp,label="naive - mddf")
plot!(N.d,N.kb_rdf,subplot=sp,label="naive - rdf")
plot!(R.d,R.kb,subplot=sp,label="new - mddf")
plot!(R.d,R.kb_rdf,subplot=sp,label="new - rdf")
plot!(legend=:topright,subplot=sp)

sp=3
plot!(ylabel="Count",subplot=sp)
plot!(N.d,N.md_count,subplot=sp,label="naive",linewidth=3,alpha=0.5)
scatter!(R.d,R.md_count,subplot=sp,label="new")
plot!(N.d,N.md_count_random,subplot=sp,label="naive - rand",linewidth=3,alpha=0.5)
scatter!(R.d,R.md_count_random,subplot=sp,label="new -rand")

sp=4
plot!(ylabel="Shell vol", subplot=sp)
plot!(N.d,N.volume.shell,subplot=sp,label="naive")
plot!(R.d,R.volume.shell,subplot=sp,label="new")

sp=5
plot!(ylabel="Sum MD", subplot=sp)
plot!(N.d,N.sum_md_count,subplot=sp,label="naive - md")
scatter!(R.d,R.sum_md_count,subplot=sp,label="new - md")

sp=6
plot!(ylabel="Sum RAND", subplot=sp)
scatter!(N.d,N.sum_md_count_random,subplot=sp,label="naive - rand")
scatter!(N.d,N.sum_rdf_count,subplot=sp,label="naive - rdf")
plot!(R.d,R.sum_md_count_random,subplot=sp,label="new - rand")
plot!(R.d,R.sum_rdf_count,subplot=sp,label="new - rdf")
plot!(legend=:topleft,subplot=sp)

sp=7
plot!(ylabel="Atom contributions", subplot=sp)
scatter!(N.d,N.solute_atom,subplot=sp,label="naive")
plot!(R.d,R.solute_atom,subplot=sp,label="new")
#scatter!(N.d,N.solvent_atom,subplot=sp,label="naive")
#plot!(R.d,R.solvent_atom,subplot=sp,label="new")

plot!(size=(800,1300))
savefig("./plots.pdf")






