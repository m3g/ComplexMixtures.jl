
using BenchmarkTools

module test
  using Dates
  import Random
  const N = 10^10

  function g1(RNG)
    return rand(RNG)
  end

  function f1()
    RNG = Random.MersenneTwister()
    s = 0.
    for i in 1:N
      s += g1(RNG)
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

println("f1: RNG inside f")
@btime test.t(test.f1)


