#
# x, y: Arrays of dimensions (nx,3), (ny,3)
# ifmol, ilmol: first and last indexes of array x to be considered (first and last atoms of molecule)
# jfmol, jlmol: first and last indexes of array y to be considered
# returns the the minimum distance between the points considered and
# the indexes of these points in x and y vectors
#

# If x is only a vector 
"""

Returns the the minimum distance between the points considered.

"""
function minimumdistance(x::T, y::Vector{T}) where {T}
    jatom = 0
    dmin = +Inf
    for (j, yj) in pairs(y)
        d = distance(x, yj)
        if d < dmin
            jatom = j
            dmin = d
        end
    end
    return dmin, 1, jatom
end

# If both are vectors of vectors

function minimumdistance(x::Vector{T}, y::Vector{T}) where {T}
    iatom = 0
    jatom = 0
    dmin = +Inf
    for (i, xi) in pairs(x)
        for (j, yj) in pairs(y)
            d = distance(xi, yj)
            if d < dmin
                iatom = i
                jatom = j
                dmin = d
            end
        end
    end
    return dmin, iatom, jatom
end

# Function that returns the distance of a reference atom as well, to be used for 
# computation of the volume shell by Monte-Caro integration

function minimumdistance(x::Vector{T}, y::Vector{T}, jrefatom::Int) where {T}
    iatom = 0
    jatom = 0
    drefatom = +Inf
    dmin = +Inf
    for (i, xi) in pairs(x)
        for (j, yj) in pairs(y)
            d = distance(xi, yj)
            # Minimum distance of any solvent atom to the solute
            if d < dmin
                iatom = i
                jatom = j
                dmin = d
            end
            # Minimum distance of the reference atom to the solute
            if j == jrefatom
                if d < drefatom
                    drefatom = d
                end
            end
        end
    end
    return dmin, iatom, jatom, drefatom
end

#
# With periodic boundary conditions
#

function minimumdistance(x::Vector{T}, y::Vector{T}, sides::T) where {T}
    iatom = 0
    jatom = 0
    dmin = +Inf
    for (i, xi) in pairs(x)
        for (j, yj) in pairs(y)
            d = distance(xi, yj, sides)
            if d < dmin
                iatom = i
                jatom = j
                dmin = d
            end
        end
    end
    return dmin, iatom, jatom
end

# Function that returns the distance of a reference atom as well, to be used for 
# computation of the volume shell by Monte-Caro integration

function minimumdistance(x::Vector{T}, y::Vector{T}, jrefatom::Int, sides::T) where {T}
    iatom = 0
    jatom = 0
    drefatom = +Inf
    dmin = +Inf
    for (i, xi) in pairs(x)
        for (j, yj) in pairs(y)
            d = distance(xi, yj, sides)
            # Minimum distance of any solvent atom to the solute
            if d < dmin
                iatom = i
                jatom = j
                dmin = d
            end
            # Minimum distance of the reference atom to the solute
            if j == jrefatom
                if d < drefatom
                    drefatom = d
                end
            end
        end
    end
    return dmin, iatom, jatom, drefatom
end
