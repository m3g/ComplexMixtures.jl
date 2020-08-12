#
# Replace default random number generator by other that does not 
# have mutiple-threading problems
#
import Random
using Future:randjump
const RNGS = [randjump(Random.MersenneTwister(),big(10)^20)]
init_random() = foreach(_ -> push!(RNGS, randjump(last(RNGS),big(10)^20)), 2:Threads.nthreads())
random() = rand(RNGS[Threads.threadid()])
random(arg) = rand(RNGS[Threads.threadid()],arg)

