"""

```
writexyz(x::Vector{T}, file::String) where T <: AbstractVector
```

Print test xyz file.

"""
function writexyz(x::Vector{T}, file::String) where {T<:AbstractVector}
    f = open(file, "w")
    nx = length(x)
    println(f, nx)
    println(f, "title")
    for i = 1:nx
        println(f, "H $(x[i][1]) $(x[i][2]) $(x[i][3])")
    end
    close(f)
    nothing
end
