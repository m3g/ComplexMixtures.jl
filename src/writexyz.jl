#
# Print test xyz file
#
function writexyz(x :: AbstractArray{Float64}, file :: String)
  f = open(file,"w")
  println(f,size(x,1))
  println(f,"title")
  for i in 1:size(x,1)
    println(f,"H $(x[i,1]) $(x[i,2]) $(x[i,3])")
  end
  close(f)
end


