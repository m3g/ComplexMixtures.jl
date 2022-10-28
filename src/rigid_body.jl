#
# Replace default random number generator by other that aleviates the
# mutiple-threading problems
#

import Random
using StableRNGs
using Future: randjump

function init_random(options)
    if options.seed > 0
        if options.StableRNG == true
            RNG = StableRNGs.StableRNG(options.seed)
        else
            RNG = [randjump(Random.MersenneTwister(options.seed), big(10)^20)]
            foreach(_ -> push!(RNG, randjump(last(RNG), big(10)^20)), 2:Threads.nthreads())
        end
    else
        if options.StableRNG == true
            seed = abs(rand(Int))
            RNG = StableRNGs.StableRNG(seed)
        else
            RNG = [randjump(Random.MersenneTwister(), big(10)^20)]
            foreach(_ -> push!(RNG, randjump(last(RNG), big(10)^20)), 2:Threads.nthreads())
        end
    end
    return RNG
end

random(RNG::StableRNGs.LehmerRNG) = rand(RNG)
random(RNG::StableRNGs.LehmerRNG, arg) = rand(RNG, arg)

random(RNG::Array{T}) where {T} = rand(RNG[Threads.threadid()])
random(RNG::Array{T}, arg) where {T} = rand(RNG[Threads.threadid()], arg)

"""
  centerofcoordinates(coor::AbstractVector{T}) where T

$(INTERNAL)

Computes the center of coordinates of a vector.

"""
function centerofcoordinates(coor::AbstractVector{T}) where {T}
    cm = zeros(T)
    for i in eachindex(coor)
        cm .= cm .+ coor[i]
    end
    cm = cm / n
    return cm
end


"""
  eulermat(beta, gamma, theta, deg::String)

$(INTERNAL)

This routine was added because it defines the rotation in the "human" way, an is thus used
to set the position of the fixed molecules. `deg` can only be `"degree"`, in which
case the angles with be considered in degrees. If no `deg` argument
is provided, radians are used.

That means: `beta` is a counterclockwise rotation around `x` axis.
            `gamma` is a counterclockwise rotation around `y` axis.
            `theta` is a counterclockwise rotation around `z` axis.


"""
function eulermat(beta, gamma, theta, deg::String)

    if deg != "degree"
        error("ERROR: to use radians just omit the last parameter")
    end

    beta = beta * π / 180
    gamma = gamma * π / 180
    theta = theta * π / 180

    return eulermat(beta, gamma, theta)
end

function eulermat(beta, gamma, theta)
    c1 = cos(beta)
    s1 = sin(beta)
    c2 = cos(gamma)
    s2 = sin(gamma)
    c3 = cos(theta)
    s3 = sin(theta)
    @SMatrix [    c2*c3           -c2*s3         s2
              (c1*s3+c3*s1*s2) (c1*c3-s1*s2*s3) -c2*s1
              (s1*s3-c1*c3*s2) (c1*s2*s3+c3*s1)  c1*c2 ]
end

"""
  random_move!(x_ref::AbstractVector{T}, 
               irefatom::Int,
               sides::T,
               x_new::AbstractVector{T}, RNG) where T

$(INTERNAL)

Function that generates a new random position for a molecule.

The new position is returned in `x`, a previously allocated array.

`x_solvent_random` might be a view of the array that contains all the solvent molecules.

"""
function random_move!(
    x_ref::AbstractVector{T},
    irefatom::Int,
    sides::T,
    x_new::AbstractVector{T},
    RNG,
) where {T}

    # To avoid boundary problems, the center of coordinates are generated in a 
    # much larger region, and wrapped aftwerwards
    scale = 100.0

    # Generate random coordiantes for the center of mass
    newcm = T(scale * (-sides[i] / 2 + random(RNG, Float64) * sides[i]) for i in 1:3)

    # Generate random rotation angles 
    beta = (2 * π) * random(RNG, Float64)
    gamma = (2 * π) * random(RNG, Float64)
    theta = (2 * π) * random(RNG, Float64)

    # Copy the coordinates of the molecule chosen to the random-coordinates vector
    @. x_new = x_ref

    # Take care that this molecule is not split by periodic boundary conditions, by
    # wrapping its coordinates around its reference atom
    wrap!(x_new, sides, x_ref[irefatom])

    # Move molecule to new position
    move!(x_new, newcm, beta, gamma, theta)

    return nothing
end

"""
  move!(x::AbstractVector, newcm::AbstractVector,beta, gamma, theta)

$(INTERNAL)

Translates and rotates a molecule according to the desired input center of coordinates and Euler rotations modifyies the vector x.

"""
function move!(x::AbstractVector, newcm::AbstractVector, beta, gamma, theta)
    cm = centerofcoordinates(x)
    A = eulermat(beta, gamma, theta)
    for i in eachindex(x)
        x[i] = A * (x[i] - cm) + newcm
    end
    return nothing
end


"""
  center_to_origin!(x::AbstractVector{T}, center::T) where {T<:AbstractVector}

Translates atoms of vectory in x array such that center is in the origin. (`x` is a vector of vectors).

"""
function center_to_origin!(x::AbstractVector{T}, center::T) where {T<:AbstractVector}
    for i in eachindex(x)
        x[i] = x[i] - center
    end
    return nothing
end



