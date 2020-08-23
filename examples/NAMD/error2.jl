
using BenchmarkTools

module test
  using Dates
  import Random
  const N = 10^10

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


