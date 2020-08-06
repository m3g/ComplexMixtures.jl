#
# Extract the contribution of a given atom type selection from teh 
# solute or solvent atomic contributions to the MDDF 
#
# s here is the solute or solvent selection
# atom_contributions is the R.solute_atom or R.solvent_atom arrays of the Result structure

# If a list of atom names is provided: 
function contrib(s :: Selection, atom_contributions :: Array{Float64}, namelist :: Vector{String})
  nbins = size(atom_contributions,1)
  natoms = size(atom_contributions,2)
  c = zeros(nbins)
  for i in 1:natoms
    if s.names[i] in namelist
      c += atom_contributions[:,i]
    end
  end
  return c
end

# If a list of indexes (within the a molecule) is provided
function contrib(s :: Selection, atom_contributions :: Array{Float64}, indexlist :: Vector{Int64})
  nbins = size(atom_contributions,1)
  natoms = size(atom_contributions,2)
  c = zeros(nbins)
  for i in 1:natoms
    if itype(i,s.natomspermol) in indexlist
      c += atom_contributions[:,i]
    end
  end
  return c
end
