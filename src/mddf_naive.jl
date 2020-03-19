#
# mddf_naive
#
# Computes the MDDF using the naive algorithm consisting of computing all distances
# useful for development purposes
#
# http://m3g.iqm.unicamp.br/
# http://github.com/m3g/MDDF
#

# Start with explicit optional input parameters

function mddf_naive(solute :: Solute,
                    solvent :: Solvent,
                    trajectory, # The type of trajectory defines the functions used to read it
                    output_name :: String,
                   ;firstframe :: Int64 = 1,
                    lastframe :: Int64 = -1,
                    stride :: Int64 = 1,
                    periodic :: Bool = true, 
                    binstep :: Float64 = 0.2,
                    irefatom :: Int64 = 1,
                    dbulk :: Float64 = 10.,
                    nintegral :: Int64 = 10,
                    cutoff :: Float64 = -1.,
                    n_random_samples :: Int64 = -1,
                    print_files :: Bool = true,
                    print_results :: Bool = true
                   )

  input = InputDetails(firstframe = firstframe,
                       lastframe = lastrame,
                       stride = stride,
                       periodic = periodic,
                       binstep = binstep,
                       irefatom = irefatom,
                       dbulk = dbulk,
                       nintegral = nintegral,
                       cutoff = cutoff,
                       n_random_samples = n_random_samples,
                       print_files = print_files,
                       print_results = print_results,
                      )

  return mddf_naive(solute,solvent,trajectory,output_name,input)

end

# start providing the InputDetail data structure

function mddf_naive(solute :: Solute,
                    solvent :: Solvent,
                    trajectory,
                    output_name :: String,
                    input :: InputDetails)

  
  # compute ibulk from dbulk (distance from which the solvent is considered bulk,
  # in the estimation of bulk density)

  if cutoff > 0. 
    usecutoff = true
    ibulk = round(Int64,dbulk/binstep) + 1
  else
    usecutoff = false
    ibulk = round(Int64,dbulk/binstep) + 1
    cutoff = dbulk
  end

  # Check if this is a single-solute or homogeneous solution situation: 
  if solute.nmols == 1
    single_solute = true
  else
    single_solute = false
  end

  # The number of random samples for numerical normalization
  if single_solute
    nsamples = n_random_samples
  else
    nsamples = round(Int64,n_random_samples/solute.nmols)
  end

  # Initializing the structure that carries all data
  
  mddf = MDDF_Data(nbins,solute,solvent,output_name;
                    firstframe,
                    lastframe,
                    stride,
                    periodic, 
                    binstep,
                    irefatom,
                    dbulk,
                    nintegral,
                    cutoff,
                    n_random_samples,
                    print_files,
                    print_results
                  )

  # Last atom to be read from the trajectory files (actually this is relevant
  # only for the NamdDCD format, I think)
  natoms = length(solute)+length(solvent)

  # Vector to annotate the molecules that belong to the bulk solution
  jmol_inbulk = Vector{Int64}(undef,solvent.nmols)

  # Vector that will contain randomly generated solvent molecules
  x_solvent_random = Array{Float64}(undef,solvent.natomspermol,3)
  
  # Total number of bulk molecules
  nbulk = 0

  # Computing all minimum-distances

  open(trajectory)
  if trajectory.nframes < lastframe
    error(" The number of frames of the trajectory is smaller than user-defined lastframe ")
  end
  for iframe in 1:lastframe

    # reading coordinates of next frame
    nextframe!(trajectory,solute,solvent)
    if iframe < firstrame 
      continue
    end
    if iframe%stride != 0
      continue
    end

    # computing the minimum distances
    for imol in 1:solute.nmols

      ifmol = (imol-1)*solute.natomspermol + 1
      ilmol = ifmol + solute.natomspermol - 1

      # associate names for code clarity
      x_solute = trajectory.x_solute
      x_solvent = tajectory.x_solvent

      # get pbc sides in this frame
      sides = getsides(trajectory,iframe)
      mddf.volume.total = mddf.volume.total + sides[1]*sides[2]*sides[3] 

      # Check if the cutoff is not too large considering the periodic cell size
      if cutoff > sides[1]/2. || cutoff > sides[2]/2. || cutoff > sides[3]/2.
        error("in MDDF: cutoff or dbulk > periodic_dimension/2 ")
      end

      # compute center of coordinates of solute to wrap solvent coordinates around it
      solute_center = centerofcoordinates(@view(x_solute[ifmol:ilmol]))
      wrap!(x_solvent,sides,solute_center)

      n_jmol_inbulk = 0
      for jmol in 1:solvent.nmols
        # first and last atoms of this solvent molecule
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1
        # Compute minimum distance 
        dmin, iatom, jatom = minimumdistance(ifmol,ilmol,x_solute,jfmol,jlmol,x_solvent)

        # Update histograms
        ibin = trunc(Int64,sqrt(dmin)/binstep)+1
        if ibin <= nbins
          mddf.count[ibin] += 1
          mddf.solute_atom[iatom,ibin] += 1 
          mddf.solvent_atom[jatom,ibin] += 1 
        else
          # computing the number of molecules in bulk, if at infinite dilution
          if usecutoff 
            if dmin < cutoff
              n_jmol_inbulk = n_jmol_inbulk + 1
              jmol_inbulk[n_jmol_inbulk] = jmol
            end
          else
            n_jmol_inbulk = n_jmol_inbulk + 1
            jmol_inbulk[n_jmol_inbulk] = jmol
          end
        end
      end # solvent molecules 
      nbulk = nbulk + n_jmol_inbulk

      #
      # Estimating the shell volumes using Monte-Carlo integration
      #
      for i in 1:nsamples
        xrnd = -sizes/2 + rand(Float64,3)*sizes/2 + solute_center
        for iatom in 1:solute.natomspermol
          iat = ifmol + iatom - 1
          d = dsquare(xrnd,x_solute,iat)
          if d < dmin
            dmin = d
          end
        end
        ibin = trunc(Int64,sqrt(dmin)/binstep)+1
        if ibin <= nbins
          mddf.volume.shell[ibin] += 1
        end
      end

      # Computing the random-solvent distribution to compute the random 
      # minimum-distance count
      for i in 1:nsamples
        # Choose randomly one molecule from the bulk
        jmol = rand(@view(jmol_inbulk[1:n_jmol_inbulk]))
        # Generate new random coordinates (translation and rotation) for this molecule
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1
        random_move!(jfmol,jlmol,x_solvent,sizes,solute_center,x_solvent_random)
        dmin, iatom, jatom = minimumdistance(ifmol,ilmol,x_solute,1,solvent.natomspermol,x_solvent_random)
        ibin = trunc(Int64,sqrt(dmin)/binstep)+1
        if ibin <= nbins
          mddf.count_random[ibin] += 1
        end
      end

      voltar

    end # solute molecules
  end # frames
  close(trajectory)

  # Averaging
  density_fix = (bulkdensity*av_totalvolume)/nrsolvent_random
  mddf.count_random = density_fix * mddf_count_random

  framecount = round(Int64,(lastframe-firstframe+1)/stride)
  mddf.count = mddf.count / framecount
  mddf.count_random = (mddf.count_random / n_random_samples) / framecount
  mddf.solute_atom = mddf.solute_atom / framecount
  mddf.solvent_atom = mddf.solvent_atom / framecount

  # Normalizing to compute distributions
  @. mddf.mddf = mddf.count / mddf.count_random
  for ibin in 1:nbins
    for i in 1:solute.natomspermol   
      mddf.solute_atom[i,ibin] = mddf.solute_atom[i,ibin] / mddf.count_random[ibin]
    end
    for j in 1:solvent.natomspermol
      mddf.solvent_atom[j,ibin] = mddf.solute_atom[j,ibin] / mddf.count_random[ibin]
    end
  end

  mddf.density.solvent = solvent.nmols / mddf.volume.total
  mddf.density.solute = solute.nmols / mddf.volume.total

  mddf.volume.shell = mddf.volume.shell / n_random_samples / framecount 
  mddf.volume.total = mddf.volume.total / framecount
  mddf.volume.bulk = mddf.volume.total - sum(mddf.volume.shell)

  if single_solute
    mddf.density.bulk = nbulk / mddf.volume.bulk
  else
    mddf.density.bulk = mddf.density.solvent
  end

  # Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol
  mole = 6.022140857e23
  convert = mole / 1.e24
voltar

  # Open output file and writes all information of this run

  if print_files
    output_files(mddf)
  end

  if print_results
    print_results(mddf)
  end

  return mddf 

end

function results(mddf)

  println()
  println(@printf("%s %12.5f","  Solvent density in simulation box (sites/A^3): ", simdensity))
  println(@printf("%s %12.5f","  Estimated bulk solvent density (sites/A^3): ", bulkdensity))
  println()                   
  println(@printf("%s %12.5f","  Molar volume of solvent in simulation box (cc/mol): ", convert/simdensity))
  println(@printf("%s %12.5f","  Molar volume of solvent in bulk (cc/mol): ", convert/bulkdensity))
  println()                   
  println(@printf("%s %12.5f","  Density scaling factor for numerical integration: ", density_fix))

  println()
  println(@printf("%s %12.5f","  Solute partial volume (cc/mol): ", solutevolume))

end

function output_files(mddf :: MDDF_Data)


  #
  # GMD computed with minimum distance
  #
  
  output = open(output_name,"w")
  open(20,file=output(1:length(output)))
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
    for i in nbins-round(Int64,1./binstep):nbins
      bulkerror = bulkerror + gmd[i]
    end
    bulkerror = bulkerror / ( round(Int64,1./binstep)+1 )
    for i in nbins-round(Int64,1./binstep):nbins
      sdbulkerror = (bulkerror - gmd[i])^2
    end do
    sdbulkerror = sqrt(sdbulkerror/(round(Int64,1./binstep)+1))
    println()
    println(@printf("  Average and standard deviation of long range (dbulk-1.) gmd: %12.5f +/- %12.5f ",bulkerror,sdbulkerror))
    println(output,@printf("#"))
    println(output,@printf("# Average and standard deviation of long range (dbulk-1.) gmd: %12.5f +/- %12.5f ",bulkerror,sdbulkerror))
  end

  ! Output table

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

    kb[i] = convert*(1./bulkdensity)*(md_sum - md_sum_random)

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

