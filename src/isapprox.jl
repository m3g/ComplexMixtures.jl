"""

Union of types to define an `isaprox` for testing.

"""
const CMTypes = Union{Result,Density,Volume,SolSummary,Options}

"""

```
Base.isapprox(r1::T, r2::T; debug=false) where T <: CMTypes
```

Function to test if two runs offered similar results. Mostly used in the 
package testing routines


"""
function Base.isapprox(r1::T, r2::T; debug = false) where {T<:CMTypes}
    check = true
    diff_list = Symbol[]
    for field in fieldnames(T)
        if field == :files
            continue
        end
        x = getfield(r1, field)
        y = getfield(r2, field)
        if !isapprox(x, y)
            check = false
            if debug
                push!(diff_list, field)
            end
        end
    end
    if debug
        for field in diff_list
            println(" Data in $field field differ. ")
        end
    end
    return check
end
