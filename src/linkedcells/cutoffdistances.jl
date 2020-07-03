#
# cutoffdistances: This routine that returns a list of the distances 
#                  between atoms that are smaller than a specified cutoff,
#                  for a given set of coordinates.
#
# L. Martinez, Sep 23, 2014. (to Julia on June, 2020)
#
# Returns nd, the number of distances smaller than the cutoff, and modifies d_in_cutoff
#

function cutoffdistances!(x_solute :: Array{Float64},
                          x_solvent :: Array{Float64},
                          lc_solute :: LinkedCells,
                          lc_solvent :: LinkedCells,
                          box :: Box, 
                          d_in_cutoff :: CutoffDistances)

 ifmol = 1
 ilmol = size(x_solute,1)
 return cutoffdistances!(ifmol,ilmol,x_solute,x_solvent,
                         lc_solute,lc_solvent,box,d_in_cutoff)

end


function cutoffdistances!(ifmol, ilmol, x_solute :: Array{Float64},
                          x_solvent :: Array{Float64},
                          lc_solute :: LinkedCells,
                          lc_solvent :: LinkedCells,
                          box :: Box, 
                          d_in_cutoff :: CutoffDistances)

  # Reset the d_in_cutoff structure 
  reset!(d_in_cutoff)

  # Wrap coordinates relative to the origin
  wrap!(x_solute,box.sides)
  wrap!(x_solvent,box.sides)

  # Compute the minimum coordinates of an atom in this cell, used to define the first cell
  box.xmin[1] = min( minimum(x_solute[:,1]), minimum(x_solvent[:,1]) )
  box.xmin[2] = min( minimum(x_solute[:,2]), minimum(x_solvent[:,2]) )
  box.xmin[3] = min( minimum(x_solute[:,3]), minimum(x_solvent[:,3]) )

  # Initialize linked lists
  initcells!(ifmol,ilmol,x_solute,box,lc_solute)
  initcells!(x_solvent,box,lc_solvent)

  # Loop over the cells
  for i in 1:box.nc[1]
    for j in 1:box.nc[2]
      for k in 1:box.nc[3]

        # First solute atom in this cell
        icell = icell1D(box.nc,i,j,k)
        iat = lc_solute.firstatom[icell]

        # Loop over solute atoms of this cell
        while iat > 0

          # Inside box
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j,k,d_in_cutoff)
        
          # Interactions of boxes that share faces

          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j,k,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j+1,k,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j,k+1,d_in_cutoff)

          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j,k,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j-1,k,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j,k-1,d_in_cutoff)
        
          # Interactions of boxes that share axes

          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j+1,k,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j,k+1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j-1,k,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j,k-1,d_in_cutoff)

          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j+1,k+1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j+1,k-1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j-1,k+1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i,j-1,k-1,d_in_cutoff)

          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j+1,k,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j,k+1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j-1,k,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j,k-1,d_in_cutoff)
                       
          # Interactions of boxes that share vertices

          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j+1,k+1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j+1,k-1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j-1,k+1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i+1,j-1,k-1,d_in_cutoff)

          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j+1,k+1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j+1,k-1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j-1,k+1,d_in_cutoff)
          cutoffdcell!(iat,x_solute,x_solvent,lc_solvent,box,i-1,j-1,k-1,d_in_cutoff)
        
          # Go to next atom of the solute in this cell
          iat = lc_solute.nextatom[iat]

        end
      end
    end
  end

end

