#
# Function to write the final results to output files
#

function write_output_files(mddf :: MDDF_Data)


  #
  # GMD computed with minimum distance
  #
  
  output = open(output_name,"w")
  println(output,@printf("#"))
  println(output,@printf("# Output of gmd.f90: Using MINIMUM distance to solute."))
  println(output,@printf("# Input file: %s", strip(inputfile)))
  println(output,@printf("# DCD file: %s", strip(dcdfile)))
  println(output,@printf("# First frame: %7i Last frame: %7i Stride: %7i",firstframe,lastframe, stride))
  println(output,@printf("# Periodic boundary conditions: "))
  println(output,@printf("#"))
  println(output,@printf("# Density of solvent in simulation box (sites/A^3): %12.5f",simdensity))
  println(output,@printf("# Density of solvent in bulk (estimated) (sites/A^3): %12.5f",bulkdensity))
  println(output,@printf("# Molar volume of solvent in simulation (cc/mol): %12.5f",convert/simdensity))
  println(output,@printf("# Molar volume of solvent in bulk (estimated) (cc/mol): %12.5f",convert/bulkdensity))
  println(output,@printf("#"))
  println(output,@printf("# Solute partial volume estimate (cc/mol): %12.5f",solutevolume))
  println(output,@printf("#"))
  println(output,@printf("# Number of atoms solute: %i7",solute.n))
  println(output,@printf("# Number of atoms of the solvent: %i7",solvent.n))
  println(output,@printf("#"))

  if usecutoff
    bulkerror = 0
    for i in ibulk:nbins
      bulkerror = bulkerror + gmd[i]
    end
    bulkerror = bulkerror / ( nbins-ibulk+1 )
    for i in ibulk:nbins
      sdbulkerror = (bulkerror - gmd[i])^2
    end
    sdbulkerror = sqrt(sdbulkerror/(nbins-ibulk+1))
    println()
    println(@printf(" Average and standard deviation of bulk-gmd: %12.5f +/- %12.5f",bulkerror,sdbulkerror))
    println(output,@printf("#"))
    println(output,@printf("# Average and standard deviation of bulk-gmd: %12.5f +/- %12.5f",bulkerror,sdbulkerror))
  else
    bulkerror = 0.
    for i in nbins-round(Int64,1/binstep):nbins
      bulkerror = bulkerror + gmd[i]
    end
    bulkerror = bulkerror / ( round(Int64,1/binstep)+1 )
    for i in nbins-round(Int64,1/binstep):nbins
      sdbulkerror = (bulkerror - gmd[i])^2
    end
    sdbulkerror = sqrt(sdbulkerror/(round(Int64,1/binstep)+1))
    println()
    println(@printf("  Average and standard deviation of long range (dbulk-1.) gmd: %12.5f +/- %12.5f ",bulkerror,sdbulkerror))
    println(output,@printf("#"))
    println(output,@printf("# Average and standard deviation of long range (dbulk-1.) gmd: %12.5f +/- %12.5f ",bulkerror,sdbulkerror))
  end

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

  md_sum = 0.
  md_sum_random = 0.
  for i in 1:nbins

    # Simple sums

    md_sum = md_sum + md_count[i]
    md_sum_random = md_sum_random + md_count_random[i]

    # KB integrals 

    kb[i] = convert*(1/bulkdensity)*(md_sum - md_sum_random)

    line = format(shellradius(i,binstep))                    #  1-DISTANCE
    line = line*"  "*format(gmd(i))                          #  2-GMD
    line = line*"  "*format(kb(i))                           #  3-KB INT
    line = line*"  "*format(md_count(i))                     #  4-MD COUNT
    line = line*"  "*format(md_count_random(i))              #  5-COUNT RAND
    line = line*"  "*format(md_sum)                          #  6-SUM MD
    line = line*"  "*format(md_sum_random)                   #  7-SUM RAND
    line = line*"  "*format(shellvolume(i))                  #  8-SHELL VOL
    println(output,line)

  end
  close(output)

  # Writting gmd per atom contributions for the solvent

  output = open(output_atom_gmd_contrib,"w")
  println(output,"# Solvent atomic contributions to total GMD. ")
  println(output,"#")
  println(output,"# Trajectory file: ", trajectory.filename)
  println(output,"#")
  println(output,"# Atoms: ")
  for i in 1:solvent.natomspermol
    println(output,@printf("#%6i  %s  %s",i, solvent.type[i], solvent.class[i]))
  end
  println(output,"#")
  string = "#     DISTANCE     GMD TOTAL" 
  for i in 1:solvent.natomspermol
    string = string*@printf("  %i12",i)
  end
  println(output,string)
  for i in 1:nbins
    string = format(shellradius(i,binstep))
    string = string*"  "*format(gmd[i])
    for j in 1:solvent.natomspermol
      string = string*"  "*format(gmd_atom_contribution[j,i])
    end
    println(output,string)
  end
  close(output)

  # Writting gmd per atom contributions for the solute

  output = open(output_atom_gmd_contrib_solute,"w")
  println(output,"# Solute atomic contributions to total GMD. ")
  println(output,"#")
  println(output,"# Trajectory file: ", trajectory.filename)
  println(output,"#")
  println(output,"# Atoms:")
  for i in 1:solute.n
    println(output,@printf("#%6i  %6i  %a  %a",i,solute.index[i],solute.type[i],soute.class[i]))
  end
  println(output,"#")
  string = "#     DISTANCE      GMD TOTAL"
  for i in 1:solute.n
    string = string*@sprintf("  %12i",i)
  end
  println(output,string)
  do i = 1, nbins
    string = format(shellradius(i,binstep))**format(gmd(i))
    for j in 1:nsolute
      string = string*"  "*format(gmd_atom_contribution_solute[j,i])
    end
    println(output,string)
  end do
  close(output)

  # Write final messages with names of output files and their content

  println()
  println(" OUTPUT FILES: ") 
  println()
  println(" Wrote solvent atomic GMD contributions to file: ", output_atom_gmd_contrib)
  println(" Wrote solute atomic GMD contributions to file: ", output_atom_gmd_contrib_solute)
  println()
  println(" Wrote main output file: ", output_name)
  println()

end

