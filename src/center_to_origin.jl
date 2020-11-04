#
# Translates atoms of vectory in x array such that center is in the origin
#
function center_to_origin!(x :: Vector{T}, center :: T) where T <: Vf3
  for i in eachindex(x)
    x[i] = x[i] - center
  end
  nothing
end
