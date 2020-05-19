#
# Structure to contain the names of the output files
#

using Parameters

@with_kw mutable struct OutputFiles

  output :: String
  solute_atoms :: String
  solvent_atoms :: String

end

