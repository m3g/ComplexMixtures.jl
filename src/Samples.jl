#
# Simple structure to contain the number of samples of each type
# of calculation to compute final results
#

struct Samples
  md :: Real
  random :: Real
end

Samples(;md :: Real = 0, random :: Real = 0) = Samples(md,random)

