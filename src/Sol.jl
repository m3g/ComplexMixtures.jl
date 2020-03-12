#
# Structure that contains the information about the solute and solvent molecules
#

struct Sol

  natoms :: Int64 # Total number of atoms
  nmols :: Int64 # Number of molecules
  natomspermol :: Int64 # Number of atoms per molecule

  index :: Vector{Int64} # Index in the coordinates of each atom
  imol :: Vector{Int64} # index of the molecule to which each atom belongs

end

# Here we associate two names just for the sake of clarity of the examples

Solute = Sol
Solvent = Sol

Solute( indexes :: Vector{Int64}; nmols :: Int64 = 0, natomspermol :: Int64 = 0 ) = 
   Sol(indexes,nmols=nmols,natomspermol=natomspemol)

Solvent( indexes :: Vector{Int64}; nmols :: Int64 = 0, natomspermol :: Int64 = 0 ) = 
   Sol(indexes,nmols=nmols,natomspermol=natomspemol)

# Function to initialize the structures

function Sol( indexes :: Vector{Int64}; nmols :: Int64 = 0, natomspermol :: Int64 = 0) 

  if nmols == 0 && natomspermol == 0
    error("In Solute, needs to set nmols or natomspermol.")
  end

  natoms = length(indexes)
  if nmols != 0
    if natoms%nmols != 0
      error(" Number of atoms in selection must be a multiple of nmols.")
    end
    natomspermol = round(Int64,n/nmols)
  else
   if natoms%natomsspermol != 0
     error(" Number of atoms in selection must be a multiple of natomspermols.")
   end
   nmols = natoms%natomspermol
  end
   
  # Setting the vector that contains the index of the molecule of each atom

  imol = Vector{Int64}(undef,natoms)
  for j in 1:nmols
    for i in 1:natomspermol
      imol[i] = j
    end
  end

  return Sol(natoms,nmols,natomspermol,indexex,imol)

end
          








