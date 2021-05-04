"""

```
which_types(s::Selection, indexes::Vector{Int})
```

Function that returns the list of the indexes of the types of the atoms in a
selection. For example, if a selection corresponds to a solvent of water molecules:
There are three types, 1, 2, and 3, corresponding to the three atoms of the
water molecule. If the indexes provided are, for instance, 11, 12, and 13, 
corresponding to a water molecule, this function will return 1, 2 and 3.

This is used to get equivalent-atom contributions to the distribution functions.
For example, the input indexes span all water molecules, the output of this
function will be still the three indexes corresponding to the three types
of atoms that exist in a water molecule. 

It is not possible to compute the contribution of *one* individual water molecule
if the distribution function was computed for all molecules. Thus, the necessity
to identify the types of atoms involved in a selection.   

"""
function which_types(s::Selection, indexes::Vector{Int}; warning=true)
  selected_types = Vector{Int}(undef,0)
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
          warning && println("WARNING: All possible types of atoms ($ntypes) of this selection were selected.")
          return selected_types
        end
      end
    end
  end
  return selected_types
end

