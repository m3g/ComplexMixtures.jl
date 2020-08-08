#
# Structure that contains the information about the solute and solvent molecules
#

struct Selection

  natoms :: Int64 # Total number of atoms
  nmols :: Int64 # Number of molecules
  natomspermol :: Int64 # Number of atoms per molecule

  index :: Vector{Int64} # Index of each atom in the full vector of coordinates
  imol :: Vector{Int64} # index of the molecule to which each atom belongs

  names :: Vector{String} # Types of the atoms, to be used in the atom-contributions

end

# Initialize providing the file name, and calling by default PDBTools.select

function Selection( file :: String, selection :: String; 
                    nmols :: Int64 = 0, natomspermol :: Int64 = 0 )
  sel = PDBTools.readPDB(file,selection)
  return Selection( sel, nmols = nmols, natomspermol = natomspermol )
end

# If the input is a vector of PDBTools.Atom types, load the index and types

function Selection( atoms :: Vector{PDBTools.Atom}; nmols :: Int64 = 0, natomspermol :: Int64 = 0)
  indexes = [ at.index for at in atoms ]
  names = [ at.name for at in atoms ]
  return Selection( indexes, names, nmols=nmols, natomspermol=natomspermol )
end

# If no names are provided, just repeat the indexes

function Selection( indexes :: Vector{Int64}; nmols :: Int64 = 0, natomspermol :: Int64 = 0)
  names = [ "$(index[i])" for i in 1:length(index) ]
  return Selection( indexes, names, nmols=nmols, natomspermol=natomspermol )
end

# Function to initialize the structures

function Selection( indexes :: Vector{Int64}, names :: Vector{String}; 
                    nmols :: Int64 = 0, natomspermol :: Int64 = 0) 

  if nmols == 0 && natomspermol == 0
    error("Set nmols or natomspermol when defining a selection.")
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
  iat = 0
  for j in 1:nmols
    for i in 1:natomspermol
      iat = iat + 1
      imol[iat] = j
    end
  end

  return Selection(natoms,nmols,natomspermol,indexes,imol,names[1:natomspermol])

end
          
import Base.show
function Base.show( io :: IO, s :: Selection )
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








