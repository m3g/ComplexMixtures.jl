#
# Function that returns the list of the types of atoms in a selection
# (that is, the indexes of the atom in the sequence of atoms of a molecule
# even if the indexes contain the indexes of many molecules of the
# same type
#
function which_types(s :: Selection, indexes :: Vector{Int64})
  selected_types = Vector{Int64}(undef,0)
  ntypes = 0
  for i in indexes
    isel = findfirst( ind -> ind == i, s.index )
    if isel == nothing
      error(" Atom in input list is not part of solvent (or solute).")
    else
      it = itype(isel,s.natomspermol)  
      if ! (it in selected_types)
        push!(selected_types,it)
        ntypes += 1
        if ntypes == s.natomspermol 
          println(" Warning: All possible types of atoms of this selection were selected.")
          return selected_types
        end
      end
    end
  end
  return selected_types
end

