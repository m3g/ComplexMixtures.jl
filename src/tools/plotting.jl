"""
    contourf_per_residue

This function requires loading the `Plots` package.

"""
function contourf_per_residue(args...; kwargs...)
    throw(ArgumentError("""\n

        contourf_per_residue could not be run. This can have two causes:

        1. This function requires loading the `Plots` package. 
           Use `using Plots` before calling this function.
        2. The function was called with the wrong arguments. 
           Please read the documentation, by typing `? contourf_per_residue`.
    
    """))
end

