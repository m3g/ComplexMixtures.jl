#
# Replace default random number generator by other that aleviates the
# mutiple-threading problems
#

import Random
using StableRNGs
using Future: randjump

function init_random(options)
    if options.seed > 0
        if options.StableRNG == true
            RNG = StableRNGs.StableRNG(options.seed)
        else
            RNG = [randjump(Random.MersenneTwister(options.seed), big(10)^20)]
            foreach(_ -> push!(RNG, randjump(last(RNG), big(10)^20)), 2:Threads.nthreads())
        end
    else
        if options.StableRNG == true
            seed = abs(rand(Int))
            RNG = StableRNGs.StableRNG(seed)
        else
            RNG = [randjump(Random.MersenneTwister(), big(10)^20)]
            foreach(_ -> push!(RNG, randjump(last(RNG), big(10)^20)), 2:Threads.nthreads())
        end
    end
    return RNG
end

random(RNG::StableRNGs.LehmerRNG) = rand(RNG)
random(RNG::StableRNGs.LehmerRNG, arg) = rand(RNG, arg)

random(RNG::Array{T}) where {T} = rand(RNG[Threads.threadid()])
random(RNG::Array{T}, arg) where {T} = rand(RNG[Threads.threadid()], arg)
