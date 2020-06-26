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

  # Number of distances smaller than the cutoff found
  nd = 0

  # Wrap coordinates relative to the origin
  wrap!(x_solute,box.sides)
  wrap!(x_solvent,box.sides)

  # Initialize linked lists
  initcells!(x_solvent,box,lc_solute)
  initcells!(x_solvent,box,lc_solvent)

  # Loop over the boxes that contain solute atoms
  index_cell_vector = 1
  icell = lc_solute.cell[index_cell_vector]
  while icell > 0
  
    # Check if this cell has a solvent atom, if not, cycle
    jcell = findfirst(jcell -> jcell == icell, lc_solvent.cell)
    if jcell == nothing
      cycle
    end

    # 3D indexes of the current cell
    i, j, k = ijkcell(box.nc,icell) 
    
    # Now, loop over the atoms of this cell, computing the distances      
    iat = lc_solute.firstatom(icell)
    while iat > 0

      # Coordinates of this solute atom
      xat = @view(x_solute[iat,1:3])

      # Inside box
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j,k,d_in_cutoff,nd)
  
      # Interactions of boxes that share faces

      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j,k,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j+1,k,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j,k+1,d_in_cutoff,nd)

      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j,k,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j-1,k,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j,k-1,d_in_cutoff,nd)
  
      # Interactions of boxes that share axes

      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j+1,k,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j,k+1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j-1,k,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j,k-1,d_in_cutoff,nd)

      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j+1,k+1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j+1,k-1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j-1,k+1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i,j-1,k-1,d_in_cutoff,nd)

      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j+1,k,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j,k+1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j-1,k,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j,k-1,d_in_cutoff,nd)
                   
      # Interactions of boxes that share vertices

      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j+1,k+1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j+1,k-1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j-1,k+1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i+1,j-1,k-1,d_in_cutoff,nd)

      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j+1,k+1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j+1,k-1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j-1,k+1,d_in_cutoff,nd)
      cutoffdcell!(iat,xat,x_solvent,lc_solvent,box,i-1,j-1,k-1,d_in_cutoff,nd)
  
      # Go to next atom of the solute in this cell
      iat = lc_solute.nextatom(iat)
    end

    # Go to next cell containing solute atoms
    index_cell_vector = index_cell_vector + 1
    icell = lc_solute.cell[index_cell_vector]
  end

  return nd
end

