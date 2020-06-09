#
# x, y: Arrays of dimensions (nx,3), (ny,3)
# ifmol, ilmol: first and last indexes of array x to be considered (first and last atoms of molecule)
# jfmol, jlmol: first and last indexes of array y to be considered
# returns the the minimum distance between the points considered and
# the indexes of these points in x and y vectors
#
function minimumdistance(ifmol :: Int64, ilmol :: Int64, x :: Array{Float64},
                         jfmol :: Int64, jlmol :: Int64, y :: Array{Float64})
  iatom = 0
  jatom = 0
  dmin = +Inf
  for i in ifmol:ilmol
     for j in jfmol:jlmol
       d = dsquare(x,y,i,j)
       if d < dmin
         iatom = i
         jatom = j
         dmin = d
       end
     end
  end
  dmin = sqrt(dmin)
  return dmin, iatom, jatom
end

# Function that returns the distance of a reference atom as well, to be used for 
# computation of the volume shell by Monte-Caro integration

function minimumdistance(ifmol :: Int64, ilmol :: Int64, x :: Array{Float64},
                         jfmol :: Int64, jlmol :: Int64, y :: Array{Float64},
                         jrefatom :: Int64)
  iatom = 0
  jatom = 0
  drefatom = 0.
  dmin = +Inf
  for i in ifmol:ilmol
     jcount = 0
     for j in jfmol:jlmol
       d = dsquare(x,y,i,j)
       if d < dmin
         iatom = i
         jatom = j
         dmin = d
       end
       jcount = jcount + 1
       if jrefatom == jcount
         drefatom = d
       end
     end
  end
  dmin = sqrt(dmin)
  drefatom = sqrt(drefatom)
  return dmin, iatom, jatom, drefatom
end

# Function to compute the minimum distance if on the input vectors
# only one atom (one set of coordinates)

function minimumdistance(x :: Vector{Float64}, jfmol, jlmol, y :: Array{Float64})
  dmin = +Inf
  for j in jfmol, jlmol
    dmin = min(dmin,dsquare(x,y,j))
  end
  return sqrt(dmin)
end




