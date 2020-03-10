#
# gmd: A program to compute gmd[1] radial distribution functions from
#      molecular dynamics simulations in NAMD DCD format.
#
#      Important: THIS IS NOT THE CLASSICAL RADIAL DISTRIBUTION
#                 FUNCTION. It is the shape-dependent RDF used
#                 for non-spherical solutes. It will only coincide
#                 with the classical RDF for perfectly spherical
#                 solutes, for instance if single atoms are used
#                 to define the solute and the solvent. 
#                 The normalization of this distribution 
#                 function is more complicated than the normalization
#                 of the radial distribution function for spherical
#                 solutes. Here, the bulk density of the solvent
#                 is estimated by the counting at long distances, and
#                 a random distribution of solvent molecules is used
#                 to estimate the volumes corresponding to each 
#                 minimum distance count. The normalization if done
#                 by dividing the actual count of sites by the
#                 expected non-interacting count estimated from this volume
#                 and the estimated bulk density.
#
# Please cite the following reference when using this package:
#
# L. Martinez, S. Shimizu, Minimum distance distribution functions for
# the analysis of the solvation of complex solutes and solvents.
# to be published.
#
# L. Martinez, Mar 13, 2014. (first version)
# Institute of Chemistry, State University of Campinas (UNICAMP)
#
# L. Martinez, Mar 1, 2017 (with KB integrals)
# Institute of Chemistry, State University of Campinas (UNICAMP)
#
# L. Martínez, Mar 20, 2020 (to Julia)
#
# http://m3g.iqm.unicamp.br/
# http://github.com/m3g/MDDF
#

using Printf

function mddf(solute :: SoluteSolvent,
              solvent :: SoluteSolvent,
              trajfile :: String,
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
              cutoff :: Float64 = -1.
              trajtype :: Type == NamdDCD
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
  if lastframe < firstframe && lastframe /= 0
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
  
  # The number of random molecules for numerical normalization 

  nrsolvent_random = nintegral*nrsolvent
  natsolvent_random = solvent.natomspermol*nrsolvent_random

  # Auxiliary solute copy for random distribution calculation

  solute2 = copy(solute)

  # Last atom to be read from the trajectory files

  natoms = length(solute)+length(solvent)
  lastatom = max(maximum(solute),maximum(solvent))

  # Maximum number of atoms in the list


  maxatom = max(natom,nsolute+max(nsolvent,natsolvent_random))

  # Initial estimate of the maximum number of small distances counted

  maxsmalld = nrsolvent_random

  # Initialization of the main data structure for the computation of small distances

  data = DistanceData(
             Frame( Vector{Float64}(undef,3),         # sides 
                    Vector{Float64}(undef,lastatom),  # x 
                    Vector{Float64}(undef,lastatom),  # y
                    Vector{Float64}(undef,lastatom)   # z
                  ),
             Groups( solute.n, solvent.n,
                     solute.index,  # indexes of solute atoms 
                     solvent.index, # indexes of solvent atoms
                   ),
             SmalldLinkedLists(
                     Vector{Int64}(undef,natoms), # iatomfirst
                     Vector{Int64}(undef,natoms), # iatomnext
                     Vector{Int64}(undef,3), # number of boxes in each dimension
                     Array{Int64}(undef,solute.n,3), # boxes of group1
                     Array{Int64}(undef,solvent.n,3), # boxes of group2
                     Vector{Float64}(undef,3), # length of box in each dimension
                     Vector{Int64}(undef,3), # maximum number of boxes in each dimension
                     cutoff
                              )
             # next structure is mutable
             SmallDistances(
                     1, # number of small distances
                     maxsmalld, # maximum number of small distances
                     Array{Int64}(undef,maxsmalld,2), # Indexes of the atoms of each distance
                     Vector{Float64}(undef,maxsmalld), # Distance
                           )
                     )

  # Arrays containing minimum-distance counts

  mind_mol = Vector{Int64}(undef,nrsolvent_random)
  mind_atom = Vector{Int64}(undef,(natsolvent_random)

  # Allocate solvent molecule (this will be used to generate random coordinates
  # for each solvent molecule, one at a time, later)
  
  solvent_molecule = Array{Float64}(undef,solvent.natomspermol,3)
  xref = Vector{Float64}(undef,solvent.natomspermol)
  yref = Vector{Float64}(undef,solvent.natomspermol)
  zref = Vector{Float64}(undef,solvent.natomspermol)
  xrnd = Vector{Float64}(undef,solvent.natomspermol)
  yrnd = Vector{Float64}(undef,solvent.natomspermol)
  zrnd = Vector{Float64}(undef,solvent.natomspermol)

  gmd_atom_contributon = Array{Float64}(undef,solvent.natomspermol,nbins)
  md_atom_contributon = Array{Float64}(undef,solvent.natomspermol,nbins)
  
  gmd_atom_contributon_solute = Array{Float64}(undef,solvent.natomspermol,nbins)
  md_atom_contributon_solute = Array{Float64}(undef,solvent.natomspermol,nbins)

  # These will contain indexes for the atoms of the randomly generated solvent molecules,
  # which are more than the number of the atoms of the solvent in the actual
  # simulation. 

  solvent_random = Vector{Int64}(undef,natsolvent_random)
  irsolv_random = Vector{Int64}(undef,natsolvent_random)

  # Opening the trajectory file, this step must return the IO stream
  # and  the number of the last frame to be read

  trajectory, lastframe = trajectory_start(trajfile)

  # Number of frames (used for normalization of counts)
  
  frames=(lastframe-firstframe+1)/stride
  println(" Number of frames to read: ", frames)

  # Counters: they are floats to avoid overflow of integers, perhaps

  md_count = zeros(Float64,nbins)
  md_count_random = zeros(Float64,nbins)
  shellvolume = zeros(Float64,nbins)
  md_atom_contribution = zeros(Float64,solvent.natomspermol,nbins)
  md_atom_contribution_solute = zeros(Float64,natoms_solute,nbins)

  bulkdensity = 0.
  simdensity = 0.
  av_totalvolume = 0.

  # Reading trajectory file and computing the gmd function
   
  for iframe in 1:nframes
   
    sides, x_solute, x_solvent = nextframe(trajectory)

    if cutoff > sides[1]/2. || cutoff > sides[2]/2. || cutoff > sides[3]/2.
      error(" ERROR in MDDF: cutoff or dbulk > periodic_dimension/2 ")
    end

    #
    # Computing the GMD data the simulation
    #

    data.frame.x_solute = x_solute
    data.frame.x_solvent = x_solvent

    # Compute all distances that are smaller than the cutoff

    smalldistances!(data)

    #
    # Computing the gmd functions from distance data
    #

    # For each solvent molecule, get the MINIMUM distance to the solute
    
    for i in 1:solvent.nmol
      mind_mol[i] = cutoff + 1. # Minimum distance found for this solvent molecule
      imind[i,1] = 0 # Solute atom corresponding to this minimum distance
      imind[i,2] = 0 # Solvent atom corresponding to this minimum distance
    end
    for i in 1:solvent.natoms
      mind_atom[i] = cutoff + 1. # Minimum distance for this atom, to compute atom contributions
    end
    for i in 1:data.smalld.n
      # Counting for computing the whole-molecule gmd 
      i1 = data.smalld.index[i,1]
      i2 = data.smalld.index[i,2]
      isolvent = solvent.imol[i2]
      if data.smalld.d[i] < mind_mol[isolvent]
        # Updating minimum distance to this solvent molecule
        mind_mol[isolvent] = data.smalld.d[i]
        # Annotating to which solute atom this md corresponds
        imind[isolvent,1] = i1
        # Annotating to which solvent atom this md corresponds
        j = i2%solvent.natomspermol
        if j == 0 
          j = solvent.natomspermol
        end
        imind[isolvent,2] = j
      end
      # Counting for computing atom-specific gmd
      if data.smalld.d[i] < mind_atom[i2]
        mind_atom[i2] = data.smalld.d[i]
      end
    end

    # Summing up current data to the gmd histogram

    for i in 1:solvent.natoms
      irad = round(Int64,nbins*mind_mol[i]/cutoff)+1
      if irad <= nbins
        # Summing up the total minimum-distance count
        md_count[irad] = md_count[irad] + 1.
        # Summing up the solvent atomic contribution
        if imind[i,2] > 0
          md_atom_contribution[imind[i,2],irad] = md_atom_contribution[imind[i,2],irad] + 1.
        end
        # Summing up the solute atomic contribution
        if imind[i,1] > 0
          md_atom_contribution_solute[imind[i,1],irad] = md_atom_contribution_solute[imind[i,1],irad] + 1.
        end
      end
    end

    # Site count at frame, to estimate the bulk density, is performed for a
    # single solvent reference site, which is taken as atom of type 'irefatom' of the solvent

    nbulk = 0
    for i in 1:solvent.natoms
      irad = round(Int64,nbins*mind_atom[i]/cutoff)+1
      if irad <= nbins
        j = i%solvent.natoms
        if j == 0 
          j = solvent.natomspermol 
        end
        if j == irefatom
          if usecutoff 
            if irad >= ibulk # Estimate bulk density from the slab between dbulk and cutoff
              nbulk = nbulk + 1
            end
          else
            if irad > nbins # Estimate from everyting beyond dbulk
              nbulk = nbulk + 1
            end
          end
        end
      end
    end

    # Total volume of the box at this frame

    totalvolume = sides[1]*sides[2]*sides[3]
    av_totalvolume = av_totalvolume + totalvolume

    # This is the average density of the solvent in the simulation box, that will
    # be averaged at the end 

    simdensity = simdensity + nrsolvent/totalvolume

    #
    # Computing random counts
    #

    # Solute coordinates are put at the first nsolute positions of x,y,z

    for i in 1:solute.natoms
      ii = solute.index[i]
      x[i] = data.frame.x[ii]
      y[i] = data.frame.y[ii]
      z[i] = data.frame.z[ii]
      solute2.index[i] = i
    end do

    #
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
      error(" ERROR: zero volume estimated for bulk region. Either the region is ",'\n'
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

  open(20,file=output_atom_gmd_contrib)
  write(20,"(a)") "# Solvent atomic contributions to total GMD. "
  write(20,"( '#',/,&
             &'# Input file: ',a,/,& 
             &'# DCD file: ',a,/,& 
             &'# Group file: ',a,/,&
             &'# PSF file: ' )")
  write(20,"(a)") "#"
  write(20,"(a)") "# Atoms: "
  do i = 1, solvent.natomspermol
    write(20,"( '#', i6, 2(tr2,a), tr2,' mass: ',f12.5 )") i, typeat(solvent(i)), classat(solvent(i)), mass(solvent(i))
  end do
  write(20,"(a)") "#"
  write(lineformat,*) "('#',t7,'DISTANCE     GMD TOTAL',",solvent.natomspermol,"(tr2,i12) )"
  write(20,lineformat) (i,i=1,solvent.natomspermol)
  write(lineformat,*) "(",solvent.natomspermol+2,"(tr2,f12.5))"
  do i = 1, nbins
    write(20,lineformat) shellradius(i,binstep), gmd(i), (gmd_atom_contribution(j,i),j=1,solvent.natomspermol)
  end do
  close(20)

  ! Writting gmd per atom contributions for the solute

  open(20,file=output_atom_gmd_contrib_solute)
  write(20,"(a)") "# Solute atomic contributions to total GMD. "
  write(20,"( '#',/,&
             &'# Input file: ',a,/,& 
             &'# DCD file: ',a,/,& 
             &'# Group file: ',a,/,&
             &'# PSF file: ' )")
  write(20,"(a)") "#"
  write(20,"(a)") "# Atoms: "
  do i = 1, nsolute
    write(20,"( '#', i6, tr2, i6, 2(tr2,a), tr2,' mass: ',f12.5 )") i, solute(i), typeat(solute(i)), &
                                                                    classat(solute(i)), mass(solute(i))
  end do
  write(20,"(a)") "#"
  write(lineformat,*) "('#',t7,'DISTANCE     GMD TOTAL',",nsolute,"(tr2,i12) )"
  write(20,lineformat) (i,i=1,nsolute)
  write(lineformat,*) "(",nsolute+2,"(tr2,f12.5))"
  do i = 1, nbins
    write(20,lineformat) shellradius(i,binstep), gmd(i), (gmd_atom_contribution_solute(j,i),j=1,nsolute)
  end do
  close(20)

  ! Write final messages with names of output files and their content
  
  time0 = etime(tarray) - time0
  write(*,*)
  write(*,"( tr2,52('-') )")
  write(*,*)
  write(*,*) ' OUTPUT FILES: ' 
  write(*,*)
  write(*,*) ' Wrote solvent atomic GMD contributions to file: ', trim(adjustl(output_atom_gmd_contrib))
  write(*,*) ' Wrote solute atomic GMD contributions to file: ', trim(adjustl(output_atom_gmd_contrib_solute))
  write(*,*)
  write(*,*) ' Wrote main output file: ', trim(adjustl(output))
  write(*,*)
  write(*,*) ' Running time: ', time0
  write(*,*) '####################################################'
  write(*,*) 
  write(*,*) '  END: Normal termination.  '
  write(*,*) 
  write(*,*) '####################################################'
  write(*,*)        

end

