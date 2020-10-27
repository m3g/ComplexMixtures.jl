#
# Replace default random number generator by other that does not 
# have mutiple-threading problems
#
import Random
#using Future:randjump
#const RNGS = [randjump(Random.MersenneTwister(),big(10)^20)]
#init_random() = foreach(_ -> push!(RNGS, randjump(last(RNGS),big(10)^20)), 2:Threads.nthreads())
#random() = rand(RNGS[Threads.threadid()])
#random(arg) = rand(RNGS[Threads.threadid()],arg)
init_random(seed) = Random.seed!(seed)
random() = rand()
random(arg) = rand(arg)

# In Julia 1.6 it seems that the generator of rand(UnitRange) will change. If that is
# the case, the tests will fail (the data has to be regenarated).
#random(arg :: UnitRange{Int64}) = arg[1] + trunc(Int64,rand()*(arg[end]-arg[1]))
