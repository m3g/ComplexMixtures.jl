#
# Structures to contain the results of the MDDF calculation
#

struct SolSummary
  natoms::Int
  nmols::Int
  natomspermol::Int
end
SolSummary(s::Selection) = SolSummary(s.natoms,s.nmols,s.natomspermol)

