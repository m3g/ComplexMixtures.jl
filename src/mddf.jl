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
# Auxiliar dimensions:
#          memory: Controls the amount of memory used for reading dcd
#                  files. If the program ends with segmentation fault
#                  without any other printing, decrease the value of
#                  this parameter.  
#
# L. Martinez, Mar 13, 2014. (first version)
# Institute of Chemistry, State University of Campinas (UNICAMP)
#
# L. Martinez, Mar 1, 2017 (with KB integrals)
# Institute of Chemistry, State University of Campinas (UNICAMP)
#
# http://leandro.iqm.unicamp.br/mdanalysis
# http://github.com/leandromaratinez98/mdanalysis
#

function mddf(solute :: Solute,
              solvent :: Solvent,
              trajfile :: String,
              output :: String,
             ;firstframe :: Int64 = 1,
              lastframe :: Int64 = -1,
              stride :: Int64 = 1,
              periodic :: Bool = true, 
              nbins :: Int64 = 1000,
              binstep :: Float64 = 0.2,
              irefatom :: Int64 = 1,
              dbulk :: Float64 = 10.,
              nint :: Int64 = 10,
              cutoff :: Float64 = 10.,
              trajtype :: Type == NamdDCD
             )

  # Arrays that are finally output

  # Constants

  pi = 4.d0*atan(1.e0)
  mole = 6.022140857e23
  memory=15000000

  sides = Array{Float64}(memory,3)
  xdcd = Array{Float64}(memory)
  ydcd = Array{Float64}(memory)
  zdcd = Array{Float64}(memory)

  ! Names of auxiliary output files
  
  output_atom_gmd_contrib = output(1:length(remove_extension(output)))//"-GMD_ATOM_CONTRIB."//file_extension(output)
  output_atom_gmd_contrib_solute = output(1:length(remove_extension(output)))//"-GMD_ATOM_SOLUTE_CONTRIB."//file_extension(output)

  ! Reading the header of psf file
  
  call getdim(psffile,inputfile,natom)
  allocate( eps(natom), sig(natom), q(natom), e(natom), s(natom), mass(natom),&
            segat(natom), resat(natom), classat(natom), typeat(natom), class(natom),&
            resid(natom), &
            irsolv(natom), solute2(natom) )

  ! Reading parameter files to get the vdW sigmas for the definition of exclusion
  ! zones
  
  nclass = 0
  open(10,file=inputfile,action='read',status='old')
  do while(.true.)
    read(10,"( a200 )",iostat=status) record
    if(status /= 0) exit
    if(keyword(record) == 'par') then
      file = value(record)
      write(*,*) ' Reading parameter file: ', file(1:length(file))
      call readpar(file,nclass,class,eps,sig)
    end if
  end do
  close(10)
  
  ! Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol

  convert = mole / 1.e24

  ! compute ibulk from dbulk (distance from which the solvent is considered bulk,
  ! in the estimation of bulk density)

  if ( dbulk-int(dbulk/binstep)*binstep > 1.e-5 ) then
    write(*,*) ' ERROR: dbulk must be a multiple of binstep. '  
    stop
  end if
  if ( usecutoff ) then
    if ( dbulk >= cutoff ) then
      write(*,*) ' ERROR: The bulk volume is zero (dbulk >= cutoff). '
      stop
    end if
    if ( (cutoff-dbulk)-int((cutoff-dbulk)/binstep)*binstep > 1.e-5 ) then
      write(*,*) ' ERROR: (cutoff-dbulk) must be a multiple of binstep. '  
      stop
    end if
    nbins = int(cutoff/binstep)
    ibulk = int(dbulk/binstep) + 1
  else
   write(*,*) ' cutoff < 0: will use ntot-n(dbulk) as nbulk '
   nbins = int(dbulk/binstep)
   ibulk = int(dbulk/binstep) + 1
   cutoff = dbulk
  end if
  write(*,*) ' Width of histogram bins: ', binstep
  write(*,*) ' Number of bins of histograms: ', nbins

  ! Allocate gmd array according to nbins

  allocate( gmd(nbins), kb(nbins), &
            md_count(nbins), md_count_random(nbins), &
            shellvolume(nbins) )
  
  ! Check for simple input errors
  
  if(stride < 1) then
    write(*,*) ' ERROR: stride cannot be less than 1. ' 
    stop
  end if
  if(lastframe < firstframe.and.lastframe /= 0) then
    write(*,*) ' ERROR: lastframe must be greater or equal to firstframe. '
    stop
  end if

  ! Output some information if not set
  
  write(*,*) ' First frame to be considered: ', firstframe
  if(lastframe == 0) then
    write(*,*) ' Last frame to be considered: last '
  else
    write(*,*) ' Last frame to be considered: ', lastframe
  end if
  write(*,*) ' Stride (will jump frames): ', stride
  write(*,*) ' Cutoff for linked cells: ', cutoff
  write(*,*) ' Bulk distance: ', dbulk
  write(*,*) ' Multiplying factor for random count: ', nintegral
  
  ! Read PSF file
  
  write(*,*) ' Reading PSF file: ', psffile(1:length(psffile))
  call readpsf(psffile,nclass,class,eps,sig,natom,segat,resat,&
               resid,classat,typeat,q,e,s,mass,.false.)
  nres = resid(natom)
  write(*,*) ' Number of atoms in PSF file: ', natom
  write(*,*) ' Number of residues in PSF file: ', nres
  write(*,*)
  
  ! Read solute and solvent information from file
  ! First reading the size of the groups to allocate arrays

  open(10,file=groupfile,action='read')
  nsolute = 0
  nsolvent = 0
  do 
    read(10,"( a200 )",iostat=status) line
    if(status /= 0) exit
    if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
      if(line(63:66) == '1.00') then      
        nsolute = nsolute + 1        
      end if
      if(line(63:66) == '2.00') then      
        nsolvent = nsolvent + 1        
      end if
    end if     
  end do
  if ( nsolute < 1 .or. nsolvent < 1 ) then
    write(*,*) ' ERROR: No atom selected for solute or solvent. '
    write(*,*) '        nsolute = ', nsolute
    write(*,*) '        nsolvent = ', nsolvent
    stop
  end if
  allocate ( solute(nsolute), solvent(nsolvent) )
  close(10)
  
  ! Now reading reading the group atoms
  
  open(10,file=groupfile,action='read')
  isolute = 0
  isolvent = 0
  mass1 = 0.
  mass2 = 0.
  iatom = 0
  do 
    read(10,"( a200 )",iostat=status) line
    if(status /= 0) exit
    if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
      iatom = iatom + 1

      ! Read atoms belonging to solute

      if(line(63:66) == '1.00') then      
        isolute = isolute + 1        
        solute(isolute) = iatom
        mass1 = mass1 + mass(iatom)
      end if
  
      ! Read atoms belonging to solvent

      if(line(63:66) == '2.00') then
        isolvent = isolvent + 1        
        solvent(isolvent) = iatom
        mass2 = mass2 + mass(iatom)
      end if

    end if     
  end do
  close(10)
  lastatom = max0(solute(nsolute),solvent(nsolvent))

  ! Counting the number of residues of the solute and solvent
  
  j = 0
  nrsolute = 0
  do i = 1, nsolute
    if(resid(solute(i)).gt.j) then
      nrsolute = nrsolute + 1 
      j = resid(solute(i))
    end if
  end do
  j = 0
  nrsolvent = 0
  do i = 1, nsolvent
    if(resid(solvent(i)).gt.j) then
      nrsolvent = nrsolvent + 1 
      j = resid(solvent(i))
    end if
    irsolv(i) = nrsolvent
  end do

  ! Output some group properties for testing purposes
  
  write(*,*) ' Number of atoms of solute: ', nsolute 
  write(*,*) ' First atom of solute: ', solute(1)
  write(*,*) ' Last atom of solute: ', solute(nsolute)
  write(*,*) ' Number of residues in solute: ', nrsolute
  write(*,*) ' Mass of solute: ', mass1
  write(*,*) ' Number of atoms of solvent: ', nsolvent 
  write(*,*) ' First atom of solvent: ', solvent(1)
  write(*,*) ' Last atom of solvent: ', solvent(nsolvent)
  write(*,*) ' Number of residues in solvent: ', nrsolvent
  write(*,*) ' Mass of solvent: ', mass2

  ! Check if the solvent atoms have obvious reading problems
  
  if ( mod(nsolvent,nrsolvent) /= 0 ) then
    write(*,*) ' ERROR: Incorrect count of solvent atoms or residues. '
    stop
  end if

  natoms_solvent = nsolvent / nrsolvent 
  write(*,*)  ' Number of atoms of each solvent molecule: ', natoms_solvent
  if ( irefatom > natoms_solvent ) then
    write(*,*) ' ERROR: Reference atom index', irefatom, ' is greater than number of '
    write(*,*) '        atoms of the solvent molecule. '
    stop
  else
    write(*,*) ' Single-site reference atom: ', irefatom,' : ', typeat(solvent(1)+irefatom-1)
  end if
  
  ! The number of random molecules for numerical normalization 

  nrsolvent_random = nintegral*nrsolvent
  natsolvent_random = natoms_solvent*nrsolvent_random

  ! Initialization of the smalldistances routine arrays
 
  maxsmalld = nrsolvent_random
  allocate( ismalld(maxsmalld,2), dsmalld(maxsmalld) )

  ! Allocate xyz and minimum-distance count arrays

  maxatom = max(natom,nsolute+max(nsolvent,natsolvent_random))
  allocate( x(maxatom), y(maxatom), z(maxatom) )
  allocate( imind(nrsolvent_random,2), mind_mol(nrsolvent_random), mind_atom(natsolvent_random) )

  ! Allocate solvent molecule (this will be used to generate random coordinates
  ! for each solvent molecule, one at a time, later)
  
  allocate( solvent_molecule(natoms_solvent,3), &
            xref(natoms_solvent), yref(natoms_solvent), zref(natoms_solvent),&
            xrnd(natoms_solvent), yrnd(natoms_solvent), zrnd(natoms_solvent) )
  allocate( gmd_atom_contribution(natoms_solvent,nbins), & 
            md_atom_contribution(natoms_solvent,nbins) )
  allocate( gmd_atom_contribution_solute(nsolute,nbins), & 
            md_atom_contribution_solute(nsolute,nbins) )

  ! These will contain indexes for the atoms of the randomly generated solvent molecules,
  ! which are more than the number of the atoms of the solvent in the actual
  ! simulation. 

  allocate( solvent_random(natsolvent_random), irsolv_random(natsolvent_random) )

  ! Checking if dcd file contains periodic cell information
  
  write(*,"( /,tr2,52('-') )")
  write(*,*) ' Periodic cell data: Read carefully. '
  call chkperiod(dcdfile,dcdaxis,readfromdcd) 
  if(.not.readfromdcd.and.periodic) then
    write(*,*) ' User provided periodic cell dimensions: '
    write(*,*) axis(1), axis(2), axis(3)
  end if
  
  ! Reading the dcd file
  
  write(*,"(tr2,52('-'),/ )")
  write(*,*) ' Reading the DCD file header: '
  open(10,file=dcdfile,action='read',form='unformatted')
  read(10) dummyc, nframes, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
  read(10) dummyi, dummyr
  read(10) ntotat
  
  write(*,*)
  write(*,*) ' Number of atoms as specified in the dcd file: ', ntotat     
  call getnframes(10,nframes,dcdaxis,lastframe)
  if( lastframe == 0 ) lastframe = nframes
  if(ntotat /= natom) then
    write(*,"(a,/,a)") ' ERROR: Number of atoms in the dcd file does not',&
                      &'        match the number of atoms in the psf file'
    stop
  end if
  
  ! Number of frames (used for normalization of counts)
  
  frames=(lastframe-firstframe+1)/stride
  write(*,*) ' Number of frames to read: ', frames

  ! Now going to read the dcd file
  
  memframes = memory / ntotat
  ncycles = (lastframe-1) / memframes + 1
  memlast = lastframe - memframes * ( ncycles - 1 )
  write(*,*) ' Will read and store in memory at most ', memframes,&
             ' frames per reading cycle. '
  write(*,*) ' There will be ', ncycles, ' cycles of reading. '
  write(*,*) ' Last cycle will read ', memlast,' frames. '
  write(*,*)        
  
  ! Reseting the counters
  
  do i = 1, nbins
    md_count(i) = 0.e0
    md_count_random(i) = 0.e0
    shellvolume(i) = 0.e0
    do j = 1, natoms_solvent
      md_atom_contribution(j,i) = 0.e0
    end do
    do j = 1, nsolute
      md_atom_contribution_solute(j,i) = 0.e0
    end do
  end do
  bulkdensity = 0.e0
  simdensity = 0.e0
  av_totalvolume = 0.e0

  ! Reading dcd file and computing the gmd function
   
  iframe = 0
  do icycle = 1, ncycles 
   
    if ( onscreenprogress ) then
      write(*,"( t3,'Cycle',t10,i5,tr2,' Reading: ',f6.2,'%')",&
            advance='no') icycle, 0. 
    end if
  
    ! Each cycle fills the memory as specified by the memory parameter 
  
    if(icycle == ncycles) then
      nfrcycle = memlast
    else
      nfrcycle = memframes
    end if
  
    iatom = 0
    do kframe = 1, nfrcycle    
      if(dcdaxis) then 
        read(10) readsidesx, t, readsidesy, t, t, readsidesz
        side(kframe,1) = sngl(readsidesx)
        side(kframe,2) = sngl(readsidesy)
        side(kframe,3) = sngl(readsidesz)
      end if
      read(10) (xdcd(k), k = iatom + 1, iatom + lastatom)
      read(10) (ydcd(k), k = iatom + 1, iatom + lastatom)            
      read(10) (zdcd(k), k = iatom + 1, iatom + lastatom)           
      iatom = iatom + ntotat
      if ( onscreenprogress ) then
        write(*,"( 7a,f6.2,'%' )",advance='no')&
             (char(8),i=1,7), 100.*float(kframe)/nfrcycle
      end if
    end do
    if ( onscreenprogress ) then
      write(*,"(' Computing: ',f6.2,'%')",advance='no') 0.
    end if
  
    ! Computing the gmd function
  
    iatom = 0
    do kframe = 1, nfrcycle
      iframe = iframe + 1

      if(mod(iframe-firstframe,stride) /= 0 .or. iframe < firstframe ) then
        iatom = iatom + ntotat
        cycle
      end if

      ! Sides of the periodic cell in this frame
  
      axis(1) = side(kframe,1) 
      axis(2) = side(kframe,2) 
      axis(3) = side(kframe,3) 

      if ( cutoff > axis(1)/2.e0 .or. &
           cutoff > axis(2)/2.e0 .or. &
           cutoff > axis(3)/2.d0 ) then
        write(*,*)
        if ( usecutoff ) then
          write(*,*) " ERROR: cutoff > periodic_dimension/2 "
        else
          write(*,*) " ERROR: dbulk > periodic_dimension/2 "
        end if
        stop
      end if

      !
      ! Computing the GMD data the simulation
      !

      do i = 1, nsolute
        ii = iatom + solute(i)
        x(solute(i)) = xdcd(ii)
        y(solute(i)) = ydcd(ii)
        z(solute(i)) = zdcd(ii)
      end do
      do i = 1, nsolvent
        ii = iatom + solvent(i)
        x(solvent(i)) = xdcd(ii)
        y(solvent(i)) = ydcd(ii)
        z(solvent(i)) = zdcd(ii)
      end do

      ! Compute all distances that are smaller than the cutoff

      memerror = .true.
      do while ( memerror ) 
        memerror = .false.
        call smalldistances(nsolute,solute,nsolvent,solvent,x,y,z,cutoff,&
                            nsmalld,ismalld,dsmalld,axis,maxsmalld,memerror)
        if ( memerror ) then
          deallocate( ismalld, dsmalld )
          maxsmalld = int(1.5*nsmalld)
          allocate( ismalld(maxsmalld,2), dsmalld(maxsmalld) )
        end if
      end do

      !
      ! Computing the gmd functions from distance data
      !

      ! For each solvent residue, get the MINIMUM distance to the solute
    
      do i = 1, nrsolvent
        mind_mol(i) = cutoff + 1.e0 ! Minimum distance found for this solvent molecule
        imind(i,1) = 0 ! Solute atom corresponding to this minimum distance
        imind(i,2) = 0 ! Solvent atom corresponding to this minimum distance
      end do
      do i = 1, nsolvent
        mind_atom(i) = cutoff + 1.e0
      end do
      do i = 1, nsmalld
        ! Counting for computing the whole-molecule gmd 
        isolvent = irsolv(ismalld(i,2))
        if ( dsmalld(i) < mind_mol(isolvent) ) then
          ! Updating minimum distance to this solvent molecule
          mind_mol(isolvent) = dsmalld(i)
          ! Annotating to which solute atom this md corresponds
          imind(isolvent,1) = ismalld(i,1)
          ! Annotating to which solvent atom this md corresponds
          j = mod(ismalld(i,2),natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          imind(isolvent,2) = j
        end if
        ! Counting for computing atom-specific gmd
        if ( dsmalld(i) < mind_atom(ismalld(i,2)) ) then
          mind_atom(ismalld(i,2)) = dsmalld(i)
        end if
      end do

      ! Summing up current data to the gmd histogram

      do i = 1, nrsolvent
        irad = int(float(nbins)*mind_mol(i)/cutoff)+1
        if( irad <= nbins ) then
          ! Summing up the total minimum-distance count
          md_count(irad) = md_count(irad) + 1.e0
          ! Summing up the solvent atomic contribution
          if ( imind(i,2) > 0 ) then
            md_atom_contribution(imind(i,2),irad) = md_atom_contribution(imind(i,2),irad) + 1.e0
          end if
          ! Summing up the solute atomic contribution
          if ( imind(i,1) > 0 ) then
            md_atom_contribution_solute(imind(i,1),irad) = md_atom_contribution_solute(imind(i,1),irad) + 1.e0
          end if
        end if
      end do

      ! Site count at frame, to estimate the bulk density, is performed for a
      ! single solvent reference site, which is taken as atom of type 'irefatom' of the solvent

      if ( usecutoff ) then
        nbulk = 0
        do i = 1, nsolvent
          irad = int(float(nbins)*mind_atom(i)/cutoff)+1
          if( irad <= nbins ) then
            j = mod(i,natoms_solvent) 
            if ( j == 0 ) j = natoms_solvent
            if ( j == irefatom .and. irad >= ibulk ) nbulk = nbulk + 1
          end if
        end do
      else
        nbulk = 0
        do i = 1, nsolvent
          irad = int(float(nbins)*mind_atom(i)/cutoff)+1
          j = mod(i,natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          if ( j == irefatom .and. irad > nbins ) nbulk = nbulk + 1
        end do
      end if

      ! Total volume of the box at this frame

      totalvolume = axis(1)*axis(2)*axis(3)
      av_totalvolume = av_totalvolume + totalvolume

      ! This is the average density of the solvent in the simulation box, that will
      ! be averaged at the end 

      simdensity = simdensity + nrsolvent/totalvolume

      !
      ! Computing random counts
      !

      ! Solute coordinates are put at the first nsolute positions of x,y,z

      do i = 1, nsolute
        ii = iatom + solute(i)
        x(i) = xdcd(ii)
        y(i) = ydcd(ii)
        z(i) = zdcd(ii)
        solute2(i) = i
      end do

      !
      ! Generating random distribution of solvent molecules in box
      !

      do isolvent_random = 1, nrsolvent_random

        ! First, pick randomly a solvent molecule from the bulk (mind_mol was just computed above
        ! for the actual simulation)
    
        ii = int((nrsolvent-1)*random()) + 1
        do while( mind_mol(ii) < dbulk ) 
          ii = int((nrsolvent-1)*random()) + 1
        end do

        ! Save the coordinates of this molecule in this frame in the solvent_molecule array
    
        jj = iatom + solvent(1) + natoms_solvent*(ii-1)
        do i = 1, natoms_solvent
          solvent_molecule(i,1) = xdcd(jj+i-1)
          solvent_molecule(i,2) = ydcd(jj+i-1)
          solvent_molecule(i,3) = zdcd(jj+i-1)
        end do

        ! Put molecule in its center of coordinates for it to be the reference coordinate for the 
        ! random coordinates that will be generated
  
        cmx = 0.
        cmy = 0.
        cmz = 0.
        do i = 1, natoms_solvent
          cmx = cmx + solvent_molecule(i,1)
          cmy = cmy + solvent_molecule(i,2)
          cmz = cmz + solvent_molecule(i,3)
        end do
        cmx = cmx / float(natoms_solvent)
        cmy = cmy / float(natoms_solvent)
        cmz = cmz / float(natoms_solvent)
        do i = 1, natoms_solvent
          xref(i) = solvent_molecule(i,1) - cmx
          yref(i) = solvent_molecule(i,2) - cmy
          zref(i) = solvent_molecule(i,3) - cmz
        end do

        ! Generate a random position for this molecule
  
        cmx = -axis(1)/2. + random()*axis(1) 
        cmy = -axis(2)/2. + random()*axis(2) 
        cmz = -axis(3)/2. + random()*axis(3) 
        beta = random()*2.e0*pi
        gamma = random()*2.e0*pi
        theta = random()*2.e0*pi
        call compcart(natoms_solvent,xref,yref,zref,xrnd,yrnd,zrnd,&
                      cmx,cmy,cmz,beta,gamma,theta)

        ! Add this molecule to x, y, z arrays

        do i = 1, natoms_solvent
          ii = nsolute + (isolvent_random-1)*natoms_solvent + i
          x(ii) = xrnd(i)
          y(ii) = yrnd(i)
          z(ii) = zrnd(i)
          solvent_random(ii-nsolute) = ii
          ! Annotate to which molecule this atom pertains
          irsolv_random(ii-nsolute) = isolvent_random
        end do

      end do

      ! The solute atom was already added to the xyz array for computing volumes, so now
      ! we have only to compute the distances

      memerror = .true.
      do while ( memerror ) 
        memerror = .false.
        call smalldistances(nsolute,solute2,natsolvent_random,solvent_random,x,y,z,cutoff,&
                            nsmalld,ismalld,dsmalld,axis,maxsmalld,memerror)
        if ( memerror ) then
          deallocate( ismalld, dsmalld )
          maxsmalld = int(1.5*nsmalld)
          allocate( ismalld(maxsmalld,2), dsmalld(maxsmalld) )
        end if
      end do

      ! Counting the number or random molecules with minimum distances
      ! in each region

      do i = 1, nrsolvent_random
        mind_mol(i) = cutoff + 1.e0
        imind(i,2) = 0
      end do
      do i = 1, nsmalld
        isolvent = irsolv_random(ismalld(i,2))
        if ( dsmalld(i) < mind_mol(isolvent) ) then
          mind_mol(isolvent) = dsmalld(i)
          j = mod(ismalld(i,2),natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          imind(isolvent,2) = j
        end if
      end do

      ! So lets count the sites at each bin distance for the non-interacting distribution 
      ! rescaled for the correct density

      do i = 1, nrsolvent_random
        irad = int(float(nbins)*mind_mol(i)/cutoff)+1
        if ( irad <= nbins ) then
          md_count_random(irad) = md_count_random(irad) + 1.e0
        end if
      end do

      ! Accumulate site count for each atom

      do i = 1, natsolvent_random
        mind_atom(i) = cutoff + 1.e0
      end do
      do i = 1, nsmalld
        if ( dsmalld(i) < mind_atom(ismalld(i,2)) ) then
          mind_atom(ismalld(i,2)) = dsmalld(i)
        end if
      end do
      if ( usecutoff ) then
        nbulk_random = 0
        do i = 1, natsolvent_random
          irad = int(float(nbins)*mind_atom(i)/cutoff)+1
          if ( irad <= nbins ) then
            j = mod(i,natoms_solvent) 
            if ( j == 0 ) j = natoms_solvent

            ! The counting of single-sites at the bulk region will be used to estimate
            ! the volumes of spherical shells of radius irad

            if ( j == irefatom ) then
              shellvolume(irad) = shellvolume(irad) + 1.e0
              if ( irad >= ibulk ) nbulk_random = nbulk_random + 1
            end if

          end if
        end do
      else
        nbulk_random = 0
        do i = 1, natsolvent_random
          irad = int(float(nbins)*mind_atom(i)/cutoff)+1
          j = mod(i,natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          if ( irad <= nbins ) then
            if ( j == irefatom ) shellvolume(irad) = shellvolume(irad) + 1.e0
          else
            if ( j == irefatom ) nbulk_random = nbulk_random + 1
          end if
        end do
      end if
      if ( nbulk_random == 0 ) then
        write(*,*) 
        write(*,*) ' ERROR: zero volume estimated for bulk region. Either the region is '
        write(*,*) '        too thin, or there is a numerical error. '
        write(*,*) ' frame = ', kframe
        stop
      end if

      ! We have just counted the number of times an atom of type 'irefatom' was found
      ! at the bulk region. The minimum-distance volume of the bulk is, then...

      bulkvolume = totalvolume*(float(nbulk_random)/nrsolvent_random)

      ! These are averaged at the end for final report:

      bulkdensity_at_frame = float(nbulk)/bulkvolume
      bulkdensity = bulkdensity + bulkdensity_at_frame

      ! Write progress

      if ( onscreenprogress ) then
        write(*,"( 7a,f6.2,'%' )",advance='no') (char(8),i=1,7), 100.*float(kframe)/nfrcycle
      else
        if ( mod(iframe,max(1,(frames/1000))) == 0 ) then
          write(*,"( '  Progress: ',f6.2,'%' )") 100.*float(iframe)/frames
        end if
      end if
  
      iatom = iatom + ntotat
    end do
    write(*,*)
  end do
  close(10)

  !
  ! Averaging results on the number of frames
  !

  bulkdensity = bulkdensity / frames
  simdensity = simdensity / frames
  av_totalvolume = av_totalvolume / frames
  density_fix = (bulkdensity*av_totalvolume)/nrsolvent_random

  write(*,*)
  write(*,"(a,f12.5)") '  Solvent density in simulation box (sites/A^3): ', simdensity
  write(*,"(a,f12.5)") '  Estimated bulk solvent density (sites/A^3): ', bulkdensity
  write(*,*)
  write(*,"(a,f12.5)") '  Molar volume of solvent in simulation box (cc/mol): ', convert/simdensity
  write(*,"(a,f12.5)") '  Molar volume of solvent in bulk (cc/mol): ', convert/bulkdensity
  write(*,*)
  write(*,"(a,f12.5)") '  Density scaling factor for numerical integration: ', density_fix

  solutevolume = convert*(bulkdensity*av_totalvolume - nrsolvent)/bulkdensity
  write(*,*)
  write(*,"(a,f12.5)") '  Solute partial volume (cc/mol): ', solutevolume
   
  do i = 1, nbins

    md_count(i) = md_count(i)/frames
    md_count_random(i) = density_fix*md_count_random(i)/frames
    do j = 1, natoms_solvent
      md_atom_contribution(j,i) = md_atom_contribution(j,i)/frames
    end do
    do j = 1, nsolute
      md_atom_contribution_solute(j,i) = md_atom_contribution_solute(j,i)/frames
    end do
    shellvolume(i) = ((shellvolume(i)/nrsolvent_random)*av_totalvolume)/frames

    ! GMD distributions

    if ( md_count_random(i) > 0.e0 ) then
      gmd(i) = md_count(i)/md_count_random(i)
      do j = 1, natoms_solvent
        gmd_atom_contribution(j,i) = md_atom_contribution(j,i)/md_count_random(i)
      end do
      do j = 1, nsolute
        gmd_atom_contribution_solute(j,i) = md_atom_contribution_solute(j,i)/md_count_random(i)
      end do
    else
      gmd(i) = 0.e0
      do j = 1, natoms_solvent
        gmd_atom_contribution(j,i) = 0.e0
      end do
      do j = 1, nsolute
        gmd_atom_contribution_solute(j,i) = 0.e0
      end do
    end if

  end do

  ! Open output file and writes all information of this run

  !
  ! GMD computed with minimum distance
  !
  
  open(20,file=output(1:length(output)))
  write(20,"( '#',/,&
             &'# Output of gmd.f90: Using MINIMUM distance to solute.',/,&
             &'# Input file: ',a,/,& 
             &'# DCD file: ',a,/,& 
             &'# Group file: ',a,/,&
             &'# PSF file: ',a,/,& 
             &'# First frame: ',i7,' Last frame: ',i7,' Stride: ',i7,/,&
             &'#',/,&
             &'# Periodic boundary conditions: ',/,&
             &'# Periodic: ',l1,' Read from DCD: ',l1,/,&
             &'#',/,&
             &'# Density of solvent in simulation box (sites/A^3): ',f12.5,/,&
             &'# Density of solvent in bulk (estimated) (sites/A^3): ',f12.5,/,&
             &'# Molar volume of solvent in simulation (cc/mol): ',f12.5,/,&
             &'# Molar volume of solvent in bulk (estimated) (cc/mol): ',f12.5,/,&
             &'#',/,&
             &'# Solute partial volume estimate (cc/mol): ',f12.5,/,&
             &'#',/,&
             &'# Number of atoms and mass of group 1: ',i7,f12.3,/,&
             &'# First and last atoms of group 1: ',i7,tr1,i7,/,&
             &'# Number of atoms and mass of group 2: ',i7,f12.3,/,&
             &'# First and last atoms of group 2: ',i7,tr1,i7,/,&
             &'#' )" )&
             &inputfile(1:length(inputfile)),&
             &dcdfile(1:length(dcdfile)),&
             &groupfile(1:length(groupfile)),&
             &psffile(1:length(psffile)),&
             &firstframe, lastframe, stride,&
             &periodic, readfromdcd, &
             &simdensity, bulkdensity, convert/simdensity, convert/bulkdensity, solutevolume, &
             &nsolute, mass1, solute(1), solute(nsolute),& 
             &nsolvent, mass2, solvent(1), solvent(nsolvent)  

  if ( usecutoff ) then
    bulkerror = 0.e0
    do i = ibulk, nbins
      bulkerror = bulkerror + gmd(i)
    end do
    bulkerror = bulkerror / ( nbins-ibulk+1 )
    do i = ibulk, nbins
      sdbulkerror = (bulkerror - gmd(i))**2
    end do
    sdbulkerror = sqrt(sdbulkerror/(nbins-ibulk+1))
    write(*,*)
    write(*,"('  Average and standard deviation of bulk-gmd: ',f12.5,' +/-',f12.5 )") bulkerror, sdbulkerror 
    write(20,"('#')")
    write(20,"('# Average and standard deviation of bulk-gmd: ',f12.5,'+/-',f12.5 )") bulkerror, sdbulkerror 
  else
    bulkerror = 0.e0
    do i = nbins-int(1./binstep), nbins
      bulkerror = bulkerror + gmd(i)
    end do
    bulkerror = bulkerror / ( int(1./binstep)+1 )
    do i = nbins-int(1./binstep), nbins
      sdbulkerror = (bulkerror - gmd(i))**2
    end do
    sdbulkerror = sqrt(sdbulkerror/(int(1./binstep)+1))
    write(*,*)
    write(*,"('  Average and standard deviation of long range (dbulk-1.) gmd: ',f12.5,' +/-',f12.5 )") bulkerror, sdbulkerror 
    write(20,"('#')")
    write(20,"('# Average and standard deviation of long range (dbulk-1.) gmd: ',f12.5,'+/-',f12.5 )") bulkerror, sdbulkerror 
  end if

  ! Output table

  write(20,"( '# COLUMNS CORRESPOND TO: ',/,&
  &'#       1  Minimum distance to solute (dmin)',/,&
  &'#       2  GMD distribution (md count normalized by md count of random-solute distribution). ',/,&
  &'#       3  Kirwood-Buff integral (cc/mol) computed [(1/bulkdensity)*(col(6)-col(7))]. ',/,&
  &'#       4  Minimum distance site count for each dmin.',/,&
  &'#       5  Minimum distance site count for each dmin for random solute distribution.',/,&
  &'#       6  Cumulative number of molecules within dmin in the simulation',/,&
  &'#       7  Cumulative number of molecules within dmin for random solute distribution.',/,&
  &'#       8  Volume of the shell of distance dmin and width binstep.')")
  write(20,"('#')") 
  write(20,"('#   1-DISTANCE         2-GMD      3-KB INT    4-MD COUNT  5-COUNT RAND      6-SUM MD    7-SUM RAND   8-SHELL VOL')")

  md_sum = 0.e0
  md_sum_random = 0.e0
  do i = 1, nbins

    ! Simple sums

    md_sum = md_sum + md_count(i)
    md_sum_random = md_sum_random + md_count_random(i)

    ! KB integrals 

    kb(i) = convert*(1.e0/bulkdensity)*(md_sum - md_sum_random)

    lineformat = "("
    lineformat = trim(adjustl(lineformat))//"tr2,"//format(shellradius(i,binstep))
    lineformat = trim(adjustl(lineformat))//",tr2,"//format(gmd(i))
    lineformat = trim(adjustl(lineformat))//",tr2,"//format(kb(i))
    lineformat = trim(adjustl(lineformat))//",tr2,"//format(md_count(i))
    lineformat = trim(adjustl(lineformat))//",tr2,"//format(md_count_random(i))
    lineformat = trim(adjustl(lineformat))//",tr2,"//format(md_sum)
    lineformat = trim(adjustl(lineformat))//",tr2,"//format(md_sum_random)
    lineformat = trim(adjustl(lineformat))//",tr2,"//format(shellvolume(i))
    lineformat = trim(adjustl(lineformat))//")"

    write(20,lineformat) &
    shellradius(i,binstep),&                                      !  1-DISTANCE
    gmd(i),&                                                      !  2-GMD
    kb(i),&                                                       !  3-KB INT
    md_count(i),&                                                 !  4-MD COUNT
    md_count_random(i),&                                          !  5-COUNT RAND
    md_sum,&                                                      !  6-SUM MD
    md_sum_random,&                                               !  7-SUM RAND
    shellvolume(i)                                                !  8-SHELL VOL

  end do
  close(20)

  ! Writting gmd per atom contributions for the solvent

  open(20,file=output_atom_gmd_contrib)
  write(20,"(a)") "# Solvent atomic contributions to total GMD. "
  write(20,"( '#',/,&
             &'# Input file: ',a,/,& 
             &'# DCD file: ',a,/,& 
             &'# Group file: ',a,/,&
             &'# PSF file: ' )")
  write(20,"(a)") "#"
  write(20,"(a)") "# Atoms: "
  do i = 1, natoms_solvent
    write(20,"( '#', i6, 2(tr2,a), tr2,' mass: ',f12.5 )") i, typeat(solvent(i)), classat(solvent(i)), mass(solvent(i))
  end do
  write(20,"(a)") "#"
  write(lineformat,*) "('#',t7,'DISTANCE     GMD TOTAL',",natoms_solvent,"(tr2,i12) )"
  write(20,lineformat) (i,i=1,natoms_solvent)
  write(lineformat,*) "(",natoms_solvent+2,"(tr2,f12.5))"
  do i = 1, nbins
    write(20,lineformat) shellradius(i,binstep), gmd(i), (gmd_atom_contribution(j,i),j=1,natoms_solvent)
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

end program g_minimum_distance

! Computes the volume of the spherical shell 
! defined within [(i-1)*step,i*step]

real function sphericalshellvolume(i,step)

  implicit none
  integer :: i
  real :: step, rmin
  real, parameter :: pi = 4.d0*atan(1.e0)
  real, parameter :: fourthirdsofpi = (4./3.)*pi

  rmin = (i-1)*step
  sphericalshellvolume = fourthirdsofpi*( (rmin+step)**3 - rmin**3 )

end function sphericalshellvolume

! Compute the point in which the radius comprises half of the
! volume of the shell

real function shellradius(i,step)

  implicit none
  integer :: i
  real :: step, rmin

  rmin = (i-1)*step
  shellradius = ( 0.5e0*( (rmin+step)**3 + rmin**3 ) )**(1.e0/3.e0)

end function shellradius

! Computes the radius that corresponds to a spherical shell of
! a given volume

real function sphereradiusfromshellvolume(volume,step)
 
  implicit none
  real :: volume, step, rmin
  real, parameter :: pi = 4.d0*atan(1.e0)
  real, parameter :: fourthirdsofpi = (4./3.)*pi
  
  if ( 3*step*volume - pi*step**4 <= 0.d0 ) then
    sphereradiusfromshellvolume = 0.d0
    return
  end if
  rmin = (sqrt(3*pi)*sqrt(3*step*volume-pi*step**4)-3*pi*step**2)/(6*pi*step)
  sphereradiusfromshellvolume = ( 0.5e0*( volume/fourthirdsofpi + 2*rmin**3 ) )**(1.e0/3.e0)

end function sphereradiusfromshellvolume

! Function to format propertly output numbers in tables

function format(x)

  implicit none
  character(len=5) :: format
  real :: x

  if ( abs(x) < 999e0 ) then
    format = 'f12.7'
  else
    format = 'e12.5'
  end if

end function format
      








