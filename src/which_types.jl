# Function that returns the list of the types of atoms in a selection
# (that is, the indexes of the atom in the sequence of atoms of a molecule
# even if the indexes contain the indexes of many molecules of the
# same type

function which_types(s :: Selection, indexes :: Vector{Int64})
  selected_types = Vector{Int64}(undef,0)
  for i in indexes
    isel = findfirst( ind -> ind == i, s.index )
    if isel == nothing
      println(" ERROR: atom in input list is not part of solvent (or solute).")
      return 0
    else
      it = itype(isel,s.natomspermol)  
      push!(selected_types,it)
    end
  end
  return selected_types
end

