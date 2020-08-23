
using BenchmarkTools

module test
  using Dates
  using Future: randjump
  import Random
  const N = 10^4

  function f0()
    s = 0.
    for i in 1:N
      s += rand()
    end
    return s
  end
  
  function f1()
    RNG1 = Random.MersenneTwister()
    s = 0.
    for i in 1:N
      s += rand(RNG1)
    end
    return s
  end
  
  const RNG = Random.MersenneTwister()
  function f2()
    s = 0.
    for i in 1:N
      s += rand(RNG)
    end
    return s
  end

  const RNGS = [randjump(Random.MersenneTwister(),big(10)^20)] 
  init_random() = foreach(_ -> push!(RNGS, randjump(last(RNGS),big(10)^20)), 2:Threads.nthreads()) 
  random() = rand(RNGS[Threads.threadid()])
  random(arg) = rand(RNGs[Threads.threadid()],arg)

  function g3()
    return random()
  end
  function f3()
    s = 0.
    for i in 1:N
      s += g3()
    end
    return s
  end

  function t(func)
    nthreads = Threads.nthreads()
    dates = Vector{DateTime}(undef,nthreads)
    task = Vector{Task}(undef,nthreads)
    for i in 1:nthreads
      task[i] = Threads.@spawn func()
      dates[i] = Dates.now()
    end
    s = 0.
    for i in 1:nthreads
      s += fetch(task[i])
    end
    #println(s)
  end

end

test.init_random()

println("f0: rand()")
@btime test.t(test.f0)

println("f1: RNG inside f")
@btime test.t(test.f1)

println("f2: using const RNG = ...")
@btime test.t(test.f2)

println("f3: using const RNGS = ...")
@btime test.t(test.f3)

println("init_random():")
@btime test.init_random()

