"""

$(TYPEDEF)

Structure to contain a list of all the distances smaller than
the cutoff.

$(TYPEDFIELDS)

"""
struct CutoffDistances

    # To store all distances smaller than the cutoff
    nd::Vector{Int} # Number of small distances (vector size = 1), just to be mutable
    d::Vector{Float64} # All distances smaller than the cutoff
    iat::Vector{Int} # atom of the solute
    jat::Vector{Int} # atom of the solvent
    imol::Vector{Int} # molecule of the solute
    jmol::Vector{Int} # molecule of the solvent

    # Size of the arrays
    maxdim::Vector{Int}

end

# Generator
CutoffDistances(natoms::Int) = CutoffDistances(
    zeros(Int, 1), # nd
    zeros(Float64, natoms), # d
    zeros(Int, natoms), # iat
    zeros(Int, natoms), # jat
    zeros(Int, natoms), # imol
    zeros(Int, natoms), # jmol
    [natoms],
) # maxdim

# Function that zeroes all the values in this structure
function reset!(c::CutoffDistances)
    c.nd[1] = 0
    @. c.d = 0.0
    @. c.iat = 0
    @. c.jat = 0
    @. c.imol = 0
    @. c.jmol = 0
    return nothing
end

function reset!(c::Vector{CutoffDistances})
    for i = 1:size(c, 1)
        reset!(c[i])
    end
    return nothing
end
