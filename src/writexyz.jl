#
# Print test xyz file
#
function writexyz(x :: AbstractArray{Float64}, file :: String)
  f = open(file,"w")
  println(f,size(x,2))
  println(f,"title")
  for i in 1:size(x,1)
    println(f,"H $(x[1,i]) $(x[2,i]) $(x[3,i])")
  end
  close(f)
  return nothing
end


