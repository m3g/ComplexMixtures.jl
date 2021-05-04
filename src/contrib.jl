"""

```
contrib(s::Selection, atom_contributions::Array{Float64}, selection)
```

Extract the contribution of a given atom type selection from the 
solute or solvent atomic contributions to the MDDF.

`s` here is the solute or solvent selection (type `ComplexMixtures.Selection`)
`atom_contributions` is the `R.solute_atom` or `R.solvent_atom` arrays of the `Result` structure,
and the last argument is the selection of atoms from the solute to be considered, given as a list 
of indexes, list of atom names, vector of `PDBTools.Atom`s, or a `PDBTools.Residue`. 

"""
function contrib(s::Selection, atom_contributions::Array{Float64}, indexes::Vector{Int})
  nbins = size(atom_contributions,1)
  c = zeros(nbins)
  for it in indexes
    if it > s.natomspermol
      error("The index list contains atoms with indexes greater than the number of atoms of the molecule.")
    end
    c += atom_contributions[:,it]
  end
  return c
end

#
# If a list of atom names is provided
#
function contrib(s::Selection, atom_contributions::Array{Float64}, names::Vector{String})
  indexes = Vector{Int}(undef,0)
  for name in names
    index = findfirst( x -> x == name, s.names )
    if index == nothing
      error(" Atom in input list is not part of solvent (or solute).")
    end
    push!(indexes,index)
  end
  return contrib(s, atom_contributions, indexes)
end

#
# If a list of atoms of PDBTools.Atom is provided
#
function contrib(s::Selection, 
                 atom_contributions::Array{Float64}, 
                 atoms::Vector{PDBTools.Atom}; 
                 warning=true)
  (warning && s.nmols > 1) && warning_nmols_types() 
  indexes = [ atom.index for atom in atoms ]
  # Check which types of atoms belong to this selection
  selected_types = which_types(s, indexes, warning=warning)
  return contrib(s, atom_contributions, selected_types)
end

#
# If a residue of type PDBTools.Residue is provided
#
function contrib(s::Selection, 
                 atom_contributions::Array{Float64}, 
                 residue::Residue;
                 warning=true)
  (warning && s.nmols > 1) && warning_nmols_types() 
  indexes = collect(residue.range)
  # Check which types of atoms belong to this selection
  selected_types = which_types(s, indexes, warning=warning)
  return contrib(s, atom_contributions, selected_types)
end

function warning_nmols_types()
  println("WARNING: there is more than one molecule in this solute. "*
          "Contributions are summed over all atoms of the same type.")
end


