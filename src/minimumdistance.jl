#
# x, y: Arrays of dimensions (nx,3), (ny,3)
# ifmol, ilmol: first and last indexes of array x to be considered (first and last atoms of molecule)
# jfmol, jlmol: first and last indexes of array y to be considered
# returns the the minimum distance between the points considered and
# the indexes of these points in x and y vectors
#

# If x is only a vector (not an array)

function minimumdistance(x :: AbstractVector{Float64}, y :: AbstractArray{Float64})
  jatom = 0
  dmin = +Inf
  ny = size(y,2)
  for j in 1:ny
    d = distance(x,@view(y[1:3,j]))
    if d < dmin
      jatom = j 
      dmin = d
    end
  end
  return dmin, 1, jatom
end

# If both are arrays

function minimumdistance(x :: AbstractArray{Float64}, y :: AbstractArray{Float64})
  iatom = 0
  jatom = 0
  dmin = +Inf
  nx = size(x,2)
  ny = size(y,2)
  for i in 1:nx
     for j in 1:ny
       d = distance(@view(x[1:3,i]),@view(y[1:3,j]))
       if d < dmin
         iatom = i 
         jatom = j 
         dmin = d
       end
     end
  end
  return dmin, iatom, jatom
end

# Function that returns the distance of a reference atom as well, to be used for 
# computation of the volume shell by Monte-Caro integration

function minimumdistance(x :: AbstractArray{Float64}, y :: AbstractArray{Float64},
                         jrefatom :: Int64)
  iatom = 0
  jatom = 0
  drefatom = +Inf
  dmin = +Inf
  nx = size(x,2)
  ny = size(y,2)
  for i in 1:nx
     for j in 1:ny
       d = distance(@view(x[1:3,i]),@view(y[1:3,j]))
       # Minimum distance of any solvent atom to the solute
       if d < dmin
         iatom = i
         jatom = j
         dmin = d
       end
       # Minimum distance of the reference atom to the solute
       if j == jrefatom
         if d < drefatom
           drefatom = d
         end
       end
     end
  end
  return dmin, iatom, jatom, drefatom
end

#
# With periodic boundary conditions
#

function minimumdistance(x :: AbstractArray{Float64}, y :: AbstractArray{Float64},
                         sides :: AbstractVector{Float64})
  iatom = 0
  jatom = 0
  dmin = +Inf
  nx = size(x,2)
  ny = size(y,2)
  for i in 1:nx
     for j in 1:ny
       d = distance(@view(x[1:3,i]),@view(y[1:3,j]),sides)
       if d < dmin
         iatom = i
         jatom = j
         dmin = d
       end
     end
  end
  return dmin, iatom, jatom
end

# Function that returns the distance of a reference atom as well, to be used for 
# computation of the volume shell by Monte-Caro integration

function minimumdistance(x :: AbstractArray{Float64}, y :: AbstractArray{Float64},
                         jrefatom :: Int64, sides :: AbstractVector{Float64})
  iatom = 0
  jatom = 0
  drefatom = +Inf
  dmin = +Inf
  nx = size(x,2)
  ny = size(y,2)
  for i in 1:nx
     for j in 1:ny
       d = distance(@view(x[1:3,i]),@view(y[1:3,j]),sides)
       # Minimum distance of any solvent atom to the solute
       if d < dmin
         iatom = i
         jatom = j
         dmin = d
       end
       # Minimum distance of the reference atom to the solute
       if j == jrefatom
         if d < drefatom
           drefatom = d
         end
       end
     end
  end
  return dmin, iatom, jatom, drefatom
end


