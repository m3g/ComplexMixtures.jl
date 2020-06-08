#
# Structure to contain the density values obtained from the calculation
#

mutable struct Density
  solute :: Float64
  solvent :: Float64
  solvent_bulk :: Float64
end

Density() = Density( 0. , 0. , 0. )

