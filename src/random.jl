#
# Replace default random number generator by other that does not 
# have mutiple-threading problems
#
import Random
random(arg) = rand(Random.MersenneTwister(),arg)
