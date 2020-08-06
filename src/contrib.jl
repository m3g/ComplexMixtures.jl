#
# Extract the contribution of a given atom type selection from teh 
# solute or solvent atomic contributions to the MDDF 
#
# s here is the solute or solvent selection
# atom_contributions is the R.solute_atom or R.solvent_atom arrays of the Result structure

#
# If a list of indexes of atoms of the simulation is provided 
#

function contrib(s :: Selection, atom_contributions :: Array{Float64}, indexes :: Vector{Int64})
  # Sum the contributions of the selected types to the contributions to be summed up    
  selected_types = which_types(s, indexes)
  nbins = size(atom_contributions,1)
  c = zeros(nbins)
  for it in selected_types
    c += atom_contributions[:,it]
  end
  return c
end

#
# Use this with, for example:
#
# atoms = PDBTools.readPDB("file.pdb")
# args = PDBTools.selindex(atoms,"resname ARG")
#

# If a list of atoms of PDBTools.Atom is provided
function contrib(s :: Selection, atom_contributions :: Array{Float64}, atoms :: Vector{PDBTools.Atom})
  indexes = [ atom.index for atom in atoms ]
  return contrib(s, atom_contributions, indexes)
end

