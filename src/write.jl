
using Statistics

function format(x)
  if abs(x) < 999
    return @sprintf("%12.7f",x)
  else
    return @sprintf("%12.5e",x)
  end
end

"""

```
write(R::Result, filename::String, solute::Selection, solvent::Selection)
```

Function to write the final results to output files as simple tables
that are human-readable and easy to analyze with other software

If the solute and solvent selections are provides, pass on the atom
names.

"""
write(R::Result, filename::String,
      solute::Selection,
      solvent::Selection) =
  write(R,filename,solute_names=solute.names,solvent_names=solvent.names)


"""

```
write(R::Result, filename::String; 
      solute_names::Vector{String} = ["nothing"], 
      solvent_names::Vector{String} = ["nothing"])
```

Optional passing of atom names.

"""
function write(R::Result, filename::String; 
               solute_names::Vector{String} = ["nothing"], 
               solvent_names::Vector{String} = ["nothing"])

  # Names of output files containing atomic contibutions
  atom_contrib_solvent = FileOperations.remove_extension(filename)*
                                 "-ATOM_CONTRIB_SOLVENT."*
                            FileOperations.file_extension(filename)
  atom_contrib_solute = FileOperations.remove_extension(filename)*
                            "-ATOM_CONTRIB_SOLUTE."*
                            FileOperations.file_extension(filename)


  #
  # GMD computed with minimum distance
  #

  # Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol
  mole = 6.022140857e23
  convert = mole / 1.e24

  output = open(filename,"w")
  println(output,@sprintf("#"))
  println(output,@sprintf("# Output of ComplexMixtures - MDDF"))
  println(output,@sprintf("# Trajectory files and weights:"))
  for i in 1:length(R.files)
    println(output,"#  $(R.files[i]) - w = $(R.weights[i])")
  end
  println(output,@sprintf("#"))
  println(output,@sprintf("# Density of solvent in simulation box (sites/A^3): %15.8f",R.density.solvent))
  println(output,@sprintf("# Density of solvent in bulk (estimated) (sites/A^3): %15.8f",R.density.solvent_bulk))
  println(output,@sprintf("# Molar volume of solvent in simulation (cc/mol): %15.8f",convert/R.density.solvent))
  println(output,@sprintf("# Molar volume of solvent in bulk (estimated) (cc/mol): %15.8f",convert/R.density.solvent_bulk))
  println(output,@sprintf("#"))
  println(output,@sprintf("# Number of atoms solute: %i9",size(R.solute_atom,2)))
  println(output,@sprintf("# Number of atoms of the solvent: %i9",size(R.solvent_atom,2)))
  println(output,@sprintf("#"))

  if R.options.usecutoff
    ibulk = setbin(R.options.dbulk,R.options.binstep)
  else
    ibulk = round(Int,R.nbins-1/R.options.binstep)
  end
  bulkerror = Statistics.mean( R.mddf[ibulk:R.nbins] )
  sdbulkerror = Statistics.std( R.mddf[ibulk:R.nbins] )
  println(output,"#")
  println(output,@sprintf("# Average and standard deviation of bulk-gmd: %12.5f +/- %12.5f",bulkerror,sdbulkerror))
  println(output,"#")

  # Output table

  println(output,"# COLUMNS CORRESPOND TO: ")
  println(output,"#       1  Minimum distance to solute (dmin)")
  println(output,"#       2  GMD distribution (md count normalized by md count of random-solute distribution)")
  println(output,"#       3  Kirwood-Buff integral (cc/mol) computed [(1/bulkdensity)*(col(6)-col(7))].")
  println(output,"#       4  Minimum distance site count for each dmin.")
  println(output,"#       5  Minimum distance site count for each dmin for random solute distribution.")
  println(output,"#       6  Cumulative number of molecules within dmin in the simulation.")
  println(output,"#       7  Cumulative number of molecules within dmin for random solute distribution.")
  println(output,"#       8  Volume of the shell of distance dmin and width binstep.")
  println(output,"#")
  println(output,"#   1-DISTANCE         2-GMD      3-KB INT    4-MD COUNT  5-COUNT RAND      6-SUM MD    7-SUM RAND   8-SHELL VOL")

  for i in 1:R.nbins
    line = "  "*format(R.d[i])                                 #  1-DISTANCE
    line = line*"  "*format(R.mddf[i])                         #  2-GMD
    line = line*"  "*format(R.kb[i])                           #  3-KB INT
    line = line*"  "*format(R.md_count[i])                     #  4-MD COUNT
    line = line*"  "*format(R.md_count_random[i])              #  5-COUNT RAND
    line = line*"  "*format(R.sum_md_count[i])                 #  6-SUM MD
    line = line*"  "*format(R.sum_md_count_random[i])          #  7-SUM RAND
    line = line*"  "*format(R.volume.shell[i])                 #  8-SHELL VOL
    println(output,line)
  end
  close(output)

  # Writting gmd per atom contributions for the solvent

  output = open(atom_contrib_solvent,"w")
  println(output,"#")
  println(output,"# Solvent atomic contributions to total MDDF. ")
  println(output,"#")
  println(output,"# Trajectory files: ", R.files)
  println(output,"#")
  println(output,"# Atoms: ")
  for i in 1:size(R.solvent_atom,2)
    if solvent_names[1] == "nothing"
      println(output,@sprintf("# %9i",i))
    else
      println(output,@sprintf("# %9i %5s",i, solvent_names[i]))
    end
  end
  println(output,"#")
  string = "#   DISTANCE     GMD TOTAL" 
  for i in 1:size(R.solvent_atom,2)
    string = string*@sprintf("  %12i",i)
  end
  println(output,string)
  for i in 1:R.nbins
    string = format(R.d[i])
    string = string*"  "*format(R.mddf[i])
    for j in 1:size(R.solvent_atom,2)
      string = string*"  "*format(R.solvent_atom[i,j])
    end
    println(output,string)
  end
  close(output)

  # Writting gmd per atom contributions for the solute

  output = open(atom_contrib_solute,"w")
  println(output,"#")
  println(output,"# Solute atomic contributions to total MDDF. ")
  println(output,"#")
  println(output,"# Trajectory files: ", R.files)
  println(output,"#")
  println(output,"# Atoms:")
  for i in 1:size(R.solute_atom,2)
    if solute_names[1] == "nothing" 
      println(output,@sprintf("# %9i",i))
    else
      println(output,@sprintf("# %9i %5s",i, solute_names[i]))
    end
  end
  println(output,"#")
  string = "#   DISTANCE      GMD TOTAL"
  for i in 1:size(R.solute_atom,2)
    string = string*@sprintf("  %12i",i)
  end
  println(output,string)
  for i in 1:R.nbins
    string = format(R.d[i])*"  "*format(R.mddf[i])
    for j in 1:size(R.solute_atom,2)
      string = string*"  "*format(R.solute_atom[i,j])
    end
    println(output,string)
  end
  close(output)

  # Write final messages with names of output files and their content

  println()
  println(" OUTPUT FILES: ") 
  println()
  println(" Wrote solvent atomic GMD contributions to file: ", atom_contrib_solvent)
  println(" Wrote solute atomic GMD contributions to file: ", atom_contrib_solute)
  println()
  println(" Wrote main output file: ", filename)
  println()

  return nothing
end

