#
# Structure to contain the density values obtained from the calculation
#

using Parameters

@with_kw mutable struct Density

  solute :: Float64 = 0.
  solvent :: Float64 = 0.
  solvent_bulk :: Float64 = 0.

end

