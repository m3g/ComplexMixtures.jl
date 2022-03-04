"""

```
center_to_origin!(x::Vector{T}, center::T) where T
```

Translates atoms of vectory in x array such that center is in the origin. (`x` is a vector of vectors).

"""
function center_to_origin!(x::Vector{T}, center::T) where {T}
    for i in eachindex(x)
        x[i] = x[i] - center
    end
    nothing
end
