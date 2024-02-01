#=

Union of types to define comparison operators.

=#
const ComplexMixturesTypes = Union{Result,Density,Volume,AtomSelection,Options}

_isapprox(x, y) = isapprox(x, y)
_isapprox(x::Vector{<:String}, y::Vector{<:String}) = all(x .== y)
_isapprox(x::String, y::String) = x == y    
_isapprox(x::Vector{<:Vector{<:Real}}, y::Vector{<:Vector{<:Real}}) = all(x .≈ y)

#=
    Base.isapprox(r1::T, r2::T; debug=false) where T <: CMTypes

Function to test if two runs offered similar results. Mostly used in the package testing routines.

=#
Base.isapprox(x::T, y::T; debug = false) where {T<:ComplexMixturesTypes} =
    _compare(_isapprox, x, y; debug = debug)

# Compare two ComplexMixtures types
import Base.==
==(x::T, y::T; debug = true) where {T<:ComplexMixturesTypes} = _compare(==, x, y; debug = debug)

function _compare(_similar::F, x::T, y::T; debug = true) where {T<:ComplexMixturesTypes} where {F<:Function}
    check = true
    diff_list = Symbol[]
    for field in fieldnames(T)
        if field in (:files, :Version)
            continue
        end
        xf = getfield(x, field)
        yf = getfield(y, field)
        try 
            if !(_similar(xf, yf))
                check = false
                if debug
                    push!(diff_list, field)
                end
            end
        catch
            error("Error comparing `$field` field of type $T.")
        end
    end
    if debug
        for field in diff_list
            println(" Data in $field field differ. ")
        end
    end
    return check
end
