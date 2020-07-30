
function f!(x :: Float64, y :: Vector{Float64})
  y[1] = 0.
  for i in 1:10^7
    y[1] = y[1] + sin(rand()*x)
  end
  return 2*y[1]
end
println("defined f")

x = 5. 
y = zeros(1)

println("x = $x, y = $y")

t = Vector{Task}(undef,2)

println("defined t = ",t)

t[1] = Threads.@spawn f!(x,y)
t[2] = Threads.@spawn f!(x,y)

check(t) = sum(@. istaskdone(t))

println("set t = ",t)

println("Before entering the loop: ", check(t))

while( check(t) < 2 )
  sleep(0.5)
  println("not done: $(check(t))")
end

println(check(t))

r = fetch(t[1]) + fetch(t[2])










