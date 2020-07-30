#
# Structure to contain the density values obtained from the calculation
#

@with_kw mutable struct Density
  solute :: Float64 = 0.
  solvent :: Float64 = 0.
  solvent_bulk :: Float64 = 0.
end

function reset!( d :: Density )
  d.solute = 0.
  d.solvent = 0.
  d.solvent_bulk = 0.
end

#function Base.show( io :: IO, d :: Density ) 
#  println(" Mean solute density: $(d.solute) ")
#  println(" Mean solvent density: $(d.solvent) ")
#  println(" Mean solvent bulk density: $(d.solvent_bulk) ")
#end


