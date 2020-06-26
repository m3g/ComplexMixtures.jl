#
# Given the index of the atom in the vector of coordinates of the solute or the solvent,
# returns the type of the atom, that is, the index of this atom within the molecule
# (goes from 1 to natomspermol)
#
function itype(iatom,natomspermol)
  itype = iatom%natomspermol
  if itype == 0
    itype = natomspermol
  end
  return itype
end
