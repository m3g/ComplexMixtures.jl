#
# smalldistances: This routine that returns a list of the distances 
#                 between atoms that are smaller than a specified cutoff,
#                 for a given set of coordinates.
#
# L. Martinez, Sep 23, 2014. (to Julia on March 3, 2020)
#
# On input: 
#    ngroup1: Number of atoms of group1 (solute)
#    group1: Vector containing the indexes in vectors x,y,z of the atoms of group1
#    ngroup2: Number of atoms of group2 (solvent) 
#    group2: Vector containing the indexes in vectors x,y,z of the atoms of group2
#    x,y,z: Coordinates of the atoms x(group1(1)), for example, is the first atom of group1
#    cutoff: The maximum distance to be considered
#    axisin: vector of dimension 3 containing the size of the periodic box in each dimension
#    maxsmalld: The estimated maximum number of distances that will be computed (see memerror)
#
# On output:
#    nsmalld: Number of distances between atoms of group1 and group1 which are smaller than 'cutoff'
#    ismalld: Indexes of the atoms of group1 and group2 corresponding to the distance reported in vector dsmalld
#             Note: the indexes go from 1 to ngroup1 and from 1 to ngroup2
#    dsmalld: The actual distance of atom of group2 to some atom of group1 (here we do not care
#             about the identity of the solute atom).
#
#    memerror: On output, it is true if the number of distances computed is greater than
#              maxsmalld. The number of distances will be stored in nsmalld, but ismalld and
#              dsmalld will not contain relevant data. If memerror is true, these two
#              vectors must be reallocated with at least "nsmalld" positions and the 
#              routine should be called again. 
#  

using OffsetArrays

# This function cycles until the vectors smalld.index and smalld.d are large
# enough to contain all the distances found

function smalldistances!( data :: DistanceData )
  memerror = true
  while memerror
    memerror = false
    smalldistances!(data,memerror)
    # If the size of arrays is not large enough, rescale them and free up the memory
    if memerror
      data.smalld.nmax = round(Int64,1.5*nsmalld)
      data.smalld.index = Array{Float64}(undef,data.smalld.nmax,2)
      data.smalld.d = Vector{Float64}(undef,data.smalld.nmax)
      GC.gc() # release memory of old arrays that were reassigned
    end
  end
end

function smalldistances!( data :: DistanceData, memerror :: Bool )

  # Associate simpler names to data variables for simplicity of the code

  nbdim = data.lists.nbdim
  iatomfirst = data.lists.iatomfirst
  iatomnext = data.lists.iatomnext
  nboxes = data.lists.nboxes
  dbox = data.lists.dbox

  sides = data.frame.sides
  x = data.frame.x
  y = data.frame.x
  z = data.frame.x

  group1 = data.groups.group1
  group2 = data.groups.group2
  group1_box = data.groups.group1_box
  group2_box = data.groups.group2_box

  cutoff = data.lists.cutoff
  cutoff2 = cutoff^2

  # Putting the atoms in their minimum image coordinates 

  center = zeros(3)
  wrap!(sides,x,y,z,center=center,sel=group1)
  wrap!(sides,x,y,z,center=center,sel=group2)

  # Prepare the linked cells
  
  xmin = -0.5 * sides
  xmax = 0.5 * sides
  
  # The side of the linked cells must be an exact divisor of the
  # box side, because of the phantom replicas 

  for i in 1:3
    boxlength = xmax[i] - xmin[i]
    nboxes[i] = max(1,boxlength/cutoff)
    dbox[i] = boxlength / nboxes[i]
  end
    
  # If the number of cells increased from the previous frame, update dimensions

  if nboxes[1] > nbdim[1] || 
     nboxes[2] > nbdim[2] || 
     nboxes[3] > nbdim[3]
    @. nbdim = nboxes
    data.lists.iatomfirst = OffsetArray{Int64}(undef,0:nboxes[1]+1,0:nboxes[2]+1,0:nboxes[3]+1)
    iatomfirst = data.lists.iatomfirst
    GC.gc()
  end

  # Reseting the iatomfirst array
  
  for i in 0:nboxes[1]+1
    for j in 0:nboxes[2]+1
      for k in 0:nboxes[3]+1
        iatomfirst[i,j,k] = 0
      end
    end
  end

  # Putting the atoms in their boxes

  for i in 1:data.ngroup1
    ibox = round(Int64,(x[group1[i]]-xmin[1])/dbox[1])+1
    jbox = round(Int64,(y[group1[i]]-xmin[2])/dbox[2])+1
    kbox = round(Int64,(z[group1[i]]-xmin[3])/dbox[3])+1
    if ibox == nbdim[1]+1 ; ibox = nbdim[1] ; end
    if jbox == nbdim[2]+1 ; jbox = nbdim[2] ; end
    if kbox == nbdim[3]+1 ; kbox = nbdim[3] ; end
    group1_box[i,1] = ibox
    group1_box[i,2] = jbox
    group1_box[i,3] = kbox
  end
  for i in 1:data.ngroup2
    ibox = round(Int64,(x[group2[i]]-xmin[1])/dbox[1])+1
    jbox = round(Int64,(y[group2[i]]-xmin[2])/dbox[2])+1
    kbox = round(Int64,(z[group2[i]]-xmin[3])/dbox[3])+1
    if ibox == nbdim[1]+1 ; ibox = nbdim[1] ; end
    if jbox == nbdim[2]+1 ; jbox = nbdim[2] ; end
    if kbox == nbdim[3]+1 ; kbox = nbdim[3] ; end
    iatomnext[i] = iatomfirst[ibox,jbox,kbox]
    iatomfirst[ibox,jbox,kbox] = i
  end

  # Filling up boundaries of periodic cell with phantom copies of the solvent atoms

  phantomcells!(nboxes,iatomfirst)

  #
  # The linked cell lists are ready, now computing the distances
  #

  # Loop over group1 atoms

  nsmalld = 0
  for igroup1 in 1:ngroup1
  
    ii = group1[igroup1]
    i = group1_box[igroup1,1] 
    j = group1_box[igroup1,2]  
    k = group1_box[igroup1,3]  

    # Check on the adjacent boxes if there is an atom of the solvent which is close enough
  
    # Inside box

    smalldcell!(data,ii,igroup1,i,j,k,memerror)
  
    # Interactions of boxes that share faces
  
    memerror = smalldcell!(data,ii,igroup1,i+1,j,k,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i,j+1,k,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i,j,k+1,memerror) 
  
    memerror = smalldcell!(data,ii,igroup1,i-1,j,k,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i,j-1,k,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i,j,k-1,memerror) 
  
    # Interactions of boxes that share axes
                 
    memerror = smalldcell!(data,ii,igroup1,i+1,j+1,k,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i+1,j,k+1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i+1,j-1,k,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i+1,j,k-1,memerror) 
  
    memerror = smalldcell!(data,ii,igroup1,i,j+1,k+1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i,j+1,k-1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i,j-1,k+1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i,j-1,k-1,memerror) 
  
    memerror = smalldcell!(data,ii,igroup1,i-1,j+1,k,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i-1,j,k+1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i-1,j-1,k,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i-1,j,k-1,memerror) 
  
    # Interactions of boxes that share vertices
  
    memerror = smalldcell!(data,ii,igroup1,i+1,j+1,k+1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i+1,j+1,k-1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i+1,j-1,k+1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i+1,j-1,k-1,memerror) 
  
    memerror = smalldcell!(data,ii,igroup1,i-1,j+1,k+1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i-1,j+1,k-1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i-1,j-1,k+1,memerror) 
    memerror = smalldcell!(data,ii,igroup1,i-1,j-1,k-1,memerror) 

  end

end
