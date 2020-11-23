# Function that sets to which histogram bin a data point pertains
# simple, but important to keep consistency over all calls

setbin(d,step) = trunc(Int, d / step ) + 1

