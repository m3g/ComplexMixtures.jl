#
# Protein - TMAO (compare new and old implementations)
#

using MDDF

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
using PDBTools
atoms = PDBTools.readPDB("./structure.pdb")

# The solute is a single protein molecule (infinte dilution case). In this case,
# use the option nmols=1
solute_indexes = [ atom.index for atom in filter( atom -> PDBTools.isprotein(atom), atoms ) ]
solute = MDDF.Selection( solute_indexes, nmols=1 )

# The solvent is TMAO, which has 14 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TMAO", atoms ) ]
solvent = MDDF.Selection( solvent_indexes, natomspermol=14 )

#solvent_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TIP3", atoms ) ]
#solvent = MDDF.Selection( solvent_indexes, natomspermol=3 )

# Alternativelly (to PDBTools, we can use VMD in background and its powerfull selections syntax,
# but you need VMD installed:
#
#solute = MDDF.Selection( MDDF.VMDselect("structure.pdb","protein",vmd="/usr/local/bin/vmd"), 
#                      nmols=1 )
#
#solvent = MDDF.Selection( MDDF.VMDselect("structure.pdb","resname TMAO",vmd="/usr/local/bin/vmd"), 
#                        natomspermol=14 ) 

# Initialize trajectroy data structure and open input stream
trajectory = MDDF.Trajectory("./trajectory.dcd",solute,solvent)

# Input options for the calcualtion
options = MDDF.Options(binstep=0.2)

# Run MDDF calculation, and get the resutls in the R structure
#R = MDDF.mddf_naive(trajectory,options)
#R = MDDF.mddf_naive(trajectory,options)
R = MDDF.mddf_linkedcells_parallel(trajectory,options)


using DelimitedFiles
using Plots
nogtk()

old = readdlm("./gmd.dat",comments=true,comment_char='#')

plot(layout=(6,1))

sp=1
plot!(ylabel="MDDF",subplot=sp)
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
savefig("./plots.pdf")
