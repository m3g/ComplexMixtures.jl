"""

$(TYPEDEF)

Structures to contain the details of a solute or solvent to
store in the results of the MDDF calculation.

$(TYPEDFIELDS)

"""
struct SolSummary
  natoms::Int
  nmols::Int
  natomspermol::Int
end
SolSummary(s::Selection) = SolSummary(s.natoms,s.nmols,s.natomspermol)

