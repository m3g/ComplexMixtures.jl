using Dates
import Random

function f()
  rng = Random.MersenneTwister()
  s = 0.
  for i in 1:(5*10^9)
    s += rand()
  end
  return s
end


function t()
  dates = Vector{DateTime}(undef,Threads.nthreads()-1)
  task = Vector{Task}(undef,Threads.nthreads()-1)
  for i in 1:Threads.nthreads()-1
    task[i] = Threads.@spawn f()
    dates[i] = Dates.now()
  end
  println(dates)
  for i in 1:Threads.nthreads()-1
    println(fetch(task[i]))
  end
end


println(Threads.nthreads())
t()

