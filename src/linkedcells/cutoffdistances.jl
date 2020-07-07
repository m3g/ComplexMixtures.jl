#
# cutoffdistances: This routine that returns a list of the distances 
#                  between atoms that are smaller than a specified cutoff,
#                  for a given set of coordinates.
#
# L. Martinez, Sep 23, 2014. (to Julia on June, 2020)
#
# Returns nd, the number of distances smaller than the cutoff, and modifies dc
#

function cutoffdistances!(x_solute :: AbstractArray{Float64},
                          x_solvent :: AbstractArray{Float64},
                          lc_solute :: LinkedCells,
                          lc_solvent :: LinkedCells,
                          box :: Box, 
                          dc :: CutoffDistances)

  # Reset the dc structure 
  reset!(dc)

  # Loop over the cells
  for i in 1:box.nc[1]
    for j in 1:box.nc[2]
      for k in 1:box.nc[3]

        # First solute atom in this cell
        icell = icell1D(box.nc,i,j,k)
        iat = lc_solute.firstatom[icell]

        # Loop over solute atoms of this cell
        while iat > 0

          # Current solute atom coordinates
          xat = @view(x_solute[iat,1:3])

          # Inside box
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j,k,dc)
        
          # Interactions of boxes that share faces

          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j,k,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j+1,k,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j,k+1,dc)

          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j,k,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j-1,k,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j,k-1,dc)
        
          # Interactions of boxes that share axes

          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j+1,k,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j,k+1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j-1,k,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j,k-1,dc)

          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j+1,k+1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j+1,k-1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j-1,k+1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j-1,k-1,dc)

          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j+1,k,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j,k+1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j-1,k,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j,k-1,dc)
                       
          # Interactions of boxes that share vertices

          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j+1,k+1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j+1,k-1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j-1,k+1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j-1,k-1,dc)

          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j+1,k+1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j+1,k-1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j-1,k+1,dc)
          cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j-1,k-1,dc)
        
          # Go to next atom of the solute in this cell
          iat = lc_solute.nextatom[iat]

        end
      end
    end
  end

end

