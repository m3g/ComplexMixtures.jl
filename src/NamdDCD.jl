# Structure to contain DCD trajectories produces with Namd
struct NamdDCD
   sides :: Vector{Float64}
   x :: Vector{Float64}
   y :: Vector{Float64}
   z :: Vector{Float64}
end
