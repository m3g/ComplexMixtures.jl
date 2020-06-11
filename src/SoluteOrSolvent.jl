#
# Structure that contains the information about the solute and solvent molecules
#

struct SoluteOrSolvent

  natoms :: Int64 # Total number of atoms
  nmols :: Int64 # Number of molecules
  natomspermol :: Int64 # Number of atoms per molecule

  index :: Vector{Int64} # Index in the coordinates of each atom
  imol :: Vector{Int64} # index of the molecule to which each atom belongs

end

# Initialize providing the file name, and calling by default VMDselect

Solute( file :: String, selection :: String; 
        vmd :: String = "vmd", nmols :: Int64 = 0, natomspermol :: Int64 = 0 ) =
  SoluteOrSolvent( VMDselect(file,selection,vmd=vmd), nmols = nmols, natomspermol = natomspermol )

Solvent( file :: String, selection :: String; 
         vmd :: String = "vmd", nmols :: Int64 = 0, natomspermol :: Int64 = 0 ) =
  SoluteOrSolvent( VMDselect(file,selection,vmd=vmd), nmols = nmols, natomspermol = natomspermol )

# Initializers from the indexes of the atoms directly

Solute( indexes :: Vector{Int64}; nmols :: Int64 = 0, natomspermol :: Int64 = 0 ) = 
   SoluteOrSolvent(indexes,nmols=nmols,natomspermol=natomspermol)

Solvent( indexes :: Vector{Int64}; nmols :: Int64 = 0, natomspermol :: Int64 = 0 ) = 
   SoluteOrSolvent(indexes,nmols=nmols,natomspermol=natomspermol)

# Function to initialize the structures

function SoluteOrSolvent( indexes :: Vector{Int64}; nmols :: Int64 = 0, natomspermol :: Int64 = 0) 

  if nmols == 0 && natomspermol == 0
    error("In Solute, needs to set nmols or natomspermol.")
  end

  natoms = length(indexes)
  if nmols != 0
    if natoms%nmols != 0
      error(" Number of atoms in selection must be a multiple of nmols.")
    end
    natomspermol = round(Int64,natoms/nmols)
  else
   if natoms%natomspermol != 0
     error(" Number of atoms in selection must be a multiple of natomspermols.")
   end
   nmols = round(Int64,natoms/natomspermol)
  end
   
  # Setting the vector that contains the index of the molecule of each atom

  imol = Vector{Int64}(undef,natoms)
  for j in 1:nmols
    for i in 1:natomspermol
      imol[i] = j
    end
  end

  return SoluteOrSolvent(natoms,nmols,natomspermol,indexes,imol)

end
          
import Base.show
function Base.show( io :: IO, s :: SoluteOrSolvent )
  if s.nmols == 1
    mol = "molecule"
  else
    mol = "molecules"
  end
  if s.natoms == 1
    at = "atom"
  else
    at = "atoms"
  end
  println(" Selection of $(s.natoms) $at belonging to $(s.nmols) $mol. ") 
end








