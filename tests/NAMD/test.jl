
module Test

  using ThreadPools

  function sumd(x)
    sumd = 0.
    for i in 1:size(x,1)
      for j in i+1:size(x,1)
        sumd += sqrt((x[i,1]-x[j,1])^2 + (x[i,2]-x[j,2])^2 + (x[i,3]-x[j,3])^2)
      end
    end
    return sumd
  end

  function serial(n)
    s = 0.
    for i in 1:Threads.nthreads()-1
      x = rand(n,3)
      s += sumd(x)
    end
    return s
  end

  function spawn(n)
    nspawn = Threads.nthreads()-1
    t = Vector{Task}(undef,nspawn)
    for i in 1:nspawn
      x = rand(n,3)
      t[i] = ThreadPools.@tspawnat i+1 sumd(x)
    end
    s = 0.
    for i in 1:nspawn
      s += fetch(t[i])
    end 
    return s
  end

end

using BenchmarkTools

n = 1000

println(" Serial: ")
@btime Test.serial(n)

println(" Parallel: ")
@btime Test.spawn(n)

