#
# Structures to contain the results of the MDDF calculation
#

struct SolSummary
  natoms :: Int64
  nmols :: Int64
  natomspermol :: Int64
end
SolSummary(s :: Selection) = SolSummary(s.natoms,s.nmols,s.natomspermol)

