"""

```
move!(x::AbstractVector{T}, newcm::AbstractVector,beta, gamma, theta) where T
```

Translates and rotates a molecule according
to the desired input center of coordinates and Euler rotations
modifyies the vector x.

"""
function move!(x::AbstractVector{T}, newcm::AbstractVector, beta, gamma, theta) where {T}

    # Compute center of coordinates of x
    cm = centerofcoordinates(x)

    # Obtain rotation matrix
    A = eulermat(beta, gamma, theta)

    for i in eachindex(x)
        x[i] = A * (x[i] - cm) + newcm
    end

    nothing
end
