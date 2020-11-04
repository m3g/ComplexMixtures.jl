#
# Print test xyz file
#
function writexyz(x :: Vector{T}, file :: String) where T <: Vf3
  f = open(file,"w")
  nx = length(x)
  println(f,nx)
  println(f,"title")
  for i in 1:nx
    println(f,"H $(x[i][1]) $(x[i][2]) $(x[i][3])")
  end
  close(f)
  nothing
end


