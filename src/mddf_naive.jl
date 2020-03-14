#
# mddf_naive
#
# Computes the MDDF using the naive algorithm consisting of computing all distances
# useful for development purposes
#
# http://m3g.iqm.unicamp.br/
# http://github.com/m3g/MDDF
#

function mddf_naive(solute :: Solute,
                    solvent :: Solvent,
                    trajectory, # The type of trajectory defines the functions used to read it
                    output_name :: String,
                   ;firstframe :: Int64 = 1,
                    lastframe :: Int64 = -1,
                    stride :: Int64 = 1,
                    periodic :: Bool = true, 
                    nbins :: Int64 = 1000,
                    binstep :: Float64 = 0.2,
                    irefatom :: Int64 = 1,
                    dbulk :: Float64 = 10.,
                    nintegral :: Int64 = 10,
                    cutoff :: Float64 = -1.,
                   )

  # Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol

  mole = 6.022140857e23
  convert = mole / 1.e24

  # Names of auxiliary output files
  output_atom_gmd_contrib = FileOperations.remove_extension(output_name)*
                            "-GMD_ATOM_CONTRIB."*
                            FileOperations.file_extension(output_name)
  output_atom_gmd_contrib_solute = FileOperations.remove_extension(output_name)*
                            "-GMD_ATOM_SOLUTE_CONTRIB."*
                            FileOperations.file_extension(output_name)

  # Check for simple input errors
  
  if stride < 1
    error(" ERROR in MDDF input: stride cannot be less than 1. ")
  end
  if lastframe < firstframe && lastframe != 0
    error(" ERROR in MDDF input: lastframe must be greater or equal to firstframe. ")
  end
  if dbulk - round(Int64,binstep*dbulk/binstep)*binstep > 1.e-5
    error("ERROR in MDDF input: dbulk must be a multiple of binstep.")
  end

  # compute ibulk from dbulk (distance from which the solvent is considered bulk,
  # in the estimation of bulk density)

  if cutoff > 0. 
    usecutoff = true
    if dbulk >= cutoff 
      error(" ERROR in MDDF input: The bulk volume is zero (dbulk must be smaller than cutoff). ")
    end
    if (cutoff-dbulk)-round(Int64,(cutoff-dbulk)/binstep)*binstep > 1.e-5 
      error(" ERROR in MDDF input: (cutoff-dbulk) must be a multiple of binstep. ")
    end
    nbins = round(Int64,cutoff/binstep)
    ibulk = round(Int64,dbulk/binstep) + 1
  else
    usecutoff = false
    println(" cutoff < 0: will use ntot-n(dbulk) as nbulk ")
    nbins = round(Int64,dbulk/binstep)
    ibulk = round(Int64,dbulk/binstep) + 1
    cutoff = dbulk
  end
  println(" Width of histogram bins: ", binstep)
  println(" Number of bins of histograms: ", nbins)

  # Allocate gmd array according to nbins

  gmd = Vector{Float64}(undef,nbins)
  kb = Vector{Float64}(undef,nbins)
  gmd_count = Vector{Int64}(undef,nbins)

  # Output some information if not set
  
  println(" First frame to be considered: ", firstframe)
  if lastframe == 0
    println(" Last frame to be considered: last ")
  else
    println(" Last frame to be considered: ", lastframe)
  end
  println(" Stride (will jump frames): ", stride)
  println(" Cutoff for linked cells: ", cutoff)
  println(" Bulk distance: ", dbulk)
  println(" Multiplying factor for random count: ", nintegral)
  
  if irefatom > solvent.natoms 
    error(" ERROR in MDDF input: Reference atom index", irefatom, " is greater than number of "*
          "                      atoms of the solvent molecule. ")
  end

  # Check if this is a single-solute or homogeneous solution situation: 
  if solute.nmols == 1
    single_solute = true
  else
    single_solute = false
  end
  
  # The number of random molecules for numerical normalization 

  nrsolvent_random = nintegral*nrsolvent
  natsolvent_random = solvent.natomspermol*nrsolvent_random

  # Auxiliary solute copy for random distribution calculation

  solute2 = copy(solute)

  # Last atom to be read from the trajectory files (actually this is relevant
  # only for the NamdDCD format, I think)

  natoms = length(solute)+length(solvent)

  # Initial estimate of the maximum number of small distances counted

  maxsmalld = nrsolvent_random

  # Initializing the structure that carries the results 
  mddf = MDDF_Results(nbins,solute,solvent)

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
      ilmol = ifmol + solute.natomspermol

      # associate names for code clarity
      x_solute = trajectory.x_solute
      x_solvent = tajectory.x_solvent

      # get pbc sides in this frame
      sides = getsides(trajectory,iframe)
      mddf.volume.total = mddf.volume.total + sides[1]*sides[2]*sides[3] 

      # Check if the cutoff is not too large considering the periodic cell size
      if cutoff > sides[1]/2. || cutoff > sides[2]/2. || cutoff > sides[3]/2.
        error(" ERROR in MDDF: cutoff or dbulk > periodic_dimension/2 ")
      end

      # compute center of coordinates of solute to wrap solvent coordinates around it
      solute_center = centerofcoordinates(@view(x_solute[ifmol:ilmol]))
      wrap!(x_solvent,sides,solute_center)

      for jmol in 1:solvent.nmols
        dmin = +Inf
        jfmol = (jmol-1)*solvent.natomspermol + 1
        mind_current_solute_atom = 0
        mind_current_solvent_atom = 0
        for iatom in 1:solute.natomspermol
           iat = ifmol + iatom - 1
           for jatom in 1:solvent.natomspermol
             jat = jfmol + jatom - 1
             d  = dsquare(x_solute,x_solvent,iat,jat)
             if d < dmin
               dmin = d
               mind_current_solute_atom = iatom
               mind_current_solvent_atom = jatom
             end

           end # solvent atoms 
        end # solute atoms

        # Update histograms
        ibin = trunc(Int64,sqrt(dmin)/binstep)+1
        if ibin <= nbins
          mddf.count[ibin] += 1
          mddf.solute_atom[mind_current_solute_atom,ibin] += 1 
          mddf.solvent_atom[mind_current_solvent_atom,ibin] += 1 
        else
          # computing the number of molecules in bulk, if at infinite dilution
          if single_solute 
            if usecutoff 
              if dmin < cutoff
                nbulk = nbulk + 1
              end
            else
              nbulk = nbulk + 1
            end
          end
        end

      end # solvent molecules 

      #
      # Estimating the shell volumes using Monte-Carlo integration
      #
      if single_solute
        for i in 1:n_random_samples
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
      end

      # Computing the random-solvent distribution to compute the random 
      # minimum-distance count

      voltar

    end # solute molecules
  end # frames
  close(trajectory)

  # Averaging
  framecount = round(Int64,(lastframe-firstframe+1)/stride)
  mddf.count = mddf.count / framecount
  mddf.count_random = (mddf.count_random / n_random_samples) / framecount
  mddf.solute_atom = mddf.solute_atom / framecount
  mddf.solvent_atom = mddf.solvent_atom / framecount

  mddf.density.solvent = solvent.nmols / mddf.volume.total
  mddf.density.solute = solute.nmols / mddf.volume.total

  mddf.volume.shell = mddf.volume.shell / n_random_samples / framecount 
  mddf.volume.total = mddf.volume.total / framecount
  mddf.volume.bulk = mddf.volume.total - sum(mddf.volume.shell)

  if single_solute
    solventdensity = solvent.nmols / volume
    bulkdensity = nbulk / bulk_volume
  else
    simdensity = solvent.nmols / volume
    bulkdensity = simdensity
  end

voltar

    # Generating random distribution of solvent molecules in box
    #

    for isolvent_random in 1:nrsolvent_random

      # First, pick randomly a solvent molecule from the bulk (mind_mol was just computed above
      # for the actual simulation)
    
      ii = round(Int64,(solvent.nmol-1)*rand(Float64))+1
      while mind_mol[ii] < dbulk 
        ii = round(Int64,(solvent.nmol-1)*rand(Float64))+1
      end

      # Save the coordinates of this molecule in this frame in the solvent_molecule array
    
      jj = solvent.index[1] + solvent.natoms*(ii-1)
      for i in 1:solvent.natoms
        solvent_molecule[i,1] = data.frame.x(jj+i-1)
        solvent_molecule[i,2] = data.frame.y(jj+i-1)
        solvent_molecule[i,3] = data.frame.z(jj+i-1)
      end

      # Put molecule in its center of coordinates for it to be the reference coordinate for the 
      # random coordinates that will be generated

      cm = centerofcoordinates(solvent_molecule)
      for i in 1:solvent.natomspermol
        xref[i] = solvent_molecule[i,1] - cm[1]
        yref[i] = solvent_molecule[i,2] - cm[2]
        zref[i] = solvent_molecule[i,3] - cm[3]
      end

      # Generate a random position for the center of coordinates of this molecule
      cm[1] = -sides[1]/2. + rand(Float64)*sides[1] 
      cm[2] = -sides[2]/2. + rand(Float64)*sides[2] 
      cm[3] = -sides[3]/2. + rand(Float64)*sides[3] 
      beta = 2*pi*rand(Float64)
      gamma = 2*pi*rand(Float64)
      theta = 2*pi*rand(Float64)
      compcart!(solvent.natomspermol,cm, 
                xref,yref,zref,beta,gamma,theta,
                xrnd, yrnd, zrnd)

      # Add this molecule to x, y, z arrays

      for i in 1:solvent.natomspermol
        ii = solute.natoms + (isolvent_random-1)*solvent.natomspermol + i
        x[ii] = xrnd[i]
        y[ii] = yrnd[i]
        z[ii] = zrnd[i]
        solvent_random[ii-nsolute] = ii
        # Annotate to which molecule this atom pertains
        irsolv_random[ii-nsolute] = isolvent_random
      end

    end

    # The solute atom was already added to the xyz array for computing volumes, so now
    # we have only to compute the distances

    smalldistances!(data)

    # Counting the number or random molecules with minimum distances
    # in each region

    for i in 1:nrsolvent_random
      mind_mol[i] = cutoff + 1.e0
      imind[i,2] = 0
    end
    for i in 1:data.smalld.n
      isolvent = irsolv_random[data.smalld.index[i,2]]
      if data.smalld.d[i] < mind_mol[isolvent]
        mind_mol[isolvent] = data.smalld.d[i]
        j = data.smalld.index[i,2]%solvent.natomspermol
        if j == 0 
          j = solvent.natomspermol
        end
        imind[isolvent,2] = j
      end
    end

    # So lets count the sites at each bin distance for the non-interacting distribution 
    # rescaled for the correct density

    for i in 1:nrsolvent_random
      irad = round(Int64,nbins*mind_mol[i]/cutoff)+1
      if irad <= nbins
        md_count_random[irad] = md_count_random[irad] + 1.
      end
    end

    # Accumulate site count for each atom

    for i in 1:natsolvent_random
      mind_atom[i] = cutoff + 1.
    end
    for i in 1:data.smalld.n
      if data.smalld.d[i] < mind_atom[data.smalld.index[i,2]]
        mind_atom[data.smalld.index[i,2]] = data.smalld.d[i]
      end
    end
    if usecutoff
      nbulk_random = 0
      for i in 1:natsolvent_random
        irad = round(Int64,nbins*mind_atom[i]/cutoff)+1
        if irad <= nbins
          j = i%solvent.natomspermol 
          if j == 0 
            j = solvent.natomspermol
          end
          # The counting of single-sites at the bulk region will be used to estimate
          # the volumes of spherical shells of radius irad
          if j == irefatom
            shellvolume[irad] = shellvolume[irad] + 1.
            if irad >= ibulk 
              nbulk_random = nbulk_random + 1
            end
          end
        end
      end
    else
      nbulk_random = 0
      for i in 1:natsolvent_random
        irad = round(Int64,nbins*mind_atom[i]/cutoff)+1
        j = i%solvent.natomspermol
        if j == 0
          j = solvent.natomspermol
        end
        if irad <= nbins
          if j == irefatom 
            shellvolume[irad] = shellvolume[irad] + 1.
          end
        else
          if j == irefatom 
            nbulk_random = nbulk_random + 1
          end
        end
      end
    end
    if nbulk_random == 0
      error(" ERROR: zero volume estimated for bulk region. Either the region is ",'\n',
            "        too thin, or there is a numerical error. ",'\n',
            " frame = ", kframe)
    end

    # We have just counted the number of times an atom of type 'irefatom' was found
    # at the bulk region. The minimum-distance volume of the bulk is, then...

    bulkvolume = totalvolume*(nbulk_random/nrsolvent_random)

    # These are averaged at the end for final report:

    bulkdensity_at_frame = nbulk/bulkvolume
    bulkdensity = bulkdensity + bulkdensity_at_frame

    iatom = iatom + ntotat
    end
    println()
  end
  trajectory_end(trajectory)

  #
  # Averaging results on the number of frames
  #

  bulkdensity = bulkdensity / frames
  simdensity = simdensity / frames
  av_totalvolume = av_totalvolume / frames
  density_fix = (bulkdensity*av_totalvolume)/nrsolvent_random

  println()
  println(@printf("%s %12.5f","  Solvent density in simulation box (sites/A^3): ", simdensity)
  println(@printf("%s %12.5f","  Estimated bulk solvent density (sites/A^3): ", bulkdensity)
  println()                   
  println(@printf("%s %12.5f","  Molar volume of solvent in simulation box (cc/mol): ", convert/simdensity)
  println(@printf("%s %12.5f","  Molar volume of solvent in bulk (cc/mol): ", convert/bulkdensity)
  println()                   
  println(@printf("%s %12.5f","  Density scaling factor for numerical integration: ", density_fix)

  solutevolume = convert*(bulkdensity*av_totalvolume - nrsolvent)/bulkdensity
  println()
  println(@printf("%s %12.5f","  Solute partial volume (cc/mol): ", solutevolume)
   
  for i in 1:nbins

    md_count[i] = md_count[i]/frames
    md_count_random[i] = density_fix*md_count_random[i]/frames
    for j in 1:solvent.natomspermol
      md_atom_contribution[j,i] = md_atom_contribution[j,i]/frames
    end
    for j in 1:solute.n
      md_atom_contribution_solute[j,i] = md_atom_contribution_solute[j,i]/frames
    end
    shellvolume[i] = ((shellvolume[i]/nrsolvent_random)*av_totalvolume)/frames

    # GMD distributions

    if md_count_random[i] > 0. 
      gmd[i] = md_count[i]/md_count_random[i]
      for j in 1:solvent.natomspermol
        gmd_atom_contribution[j,i] = md_atom_contribution[j,i]/md_count_random[i]
      end
      for j in 1:nsolute
        gmd_atom_contribution_solute[j,i] = md_atom_contribution_solute[j,i]/md_count_random[i]
      end
    else
      gmd[i] = 0.e0
      for j in 1:solvent.natomspermol
        gmd_atom_contribution[j,i] = 0.
      end
      for j in 1:nsolute
        gmd_atom_contribution_solute[j,i] = 0.
      end
    end

  end

  # Open output file and writes all information of this run

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

