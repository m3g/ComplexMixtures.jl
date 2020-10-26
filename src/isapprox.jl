#
# Function to test if two runs offered similar results. Mostly used in the 
# package testing routines
#

import Base.isapprox
const CMTypes = Union{Result,Density,Volume,SolSummary,Options} 
function isapprox( r1 :: T, r2 :: T; debug = false ) where T <: CMTypes
  check = true
  for field in fieldnames(typeof(r1)) 
    if field == :files
      continue
    end
    x = getfield(r1,field)
    y = getfield(r2,field)
    if ! isapprox(x,y)
      check = false
      if debug
        println(" Data in $field field differ. ")
      end
    end
  end
  return check
end

