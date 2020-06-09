# Function that sets to which histogram bin a data point pertains
# simple, but important to keep consistency over all calls

setbin(d :: Float64, step :: Float64) = trunc(Int64, d / step ) + 1

