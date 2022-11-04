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

@testitem "eulermat" begin
    using ComplexMixtures
    @test ComplexMixtures.eulermat(0.0, 0.0, 0.0) ≈ [1 0 0 ; 0 1 0 ; 0 0 1]
    @test ComplexMixtures.eulermat(π, 0.0, 0.0) ≈ [1 0 0 ; 0 -1 0 ; 0 0 -1]
    @test ComplexMixtures.eulermat(0.0, π, 0.0) ≈ [-1 0 0 ; 0 1 0 ; 0 0 -1]
    @test ComplexMixtures.eulermat(0.0, 0.0, π) ≈ [-1 0 0 ; 0 -1 0 ; 0 0 1]
end

"""
  center_to_origin!(x::AbstractVector{T}, center::T) where {T<:AbstractVector}

$(INTERNAL)

Translates atoms of vectory in x array such that center is in the origin. (`x` is a vector of vectors).

"""
function center_to_origin!(x::AbstractVector{T}, center::T) where {T<:SVector}
    for i in eachindex(x)
        x[i] = x[i] - center
    end
    return x
end

"""
    wrap!(x::AbstractVector{T}, vref, box::CellListMap.Box) where T<:SVector

$(INTERNAL)

Wrap the coordinates of a molecule such that its atoms are not spread in different images
of the periodic box. This is needed to generate correct random rotations for the molecule.

"""
function wrap!(x::AbstractVector{T}, vref, box::CellListMap.Box) where T<:SVector
    for i in eachindex(x)
        x[i] = CellListMap.wrap_relative_to(x[i], vref, box)
    end
    return x
end

@testitem "wrap!" begin
    using ComplexMixtures
    using StaticArrays
    import CellListMap
    # Orthorhombic box
    box = CellListMap.Box(SVector(10.0, 10.0, 10.0), 1.0)
    x = [ SVector(1.0, 0.0, 0.0), SVector(10.0, 0.0, 0.0)]
    @test ComplexMixtures.wrap!(x, x[1], box) ≈ SVector{3, Float64}[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]] 
    x = [ SVector(1.0, 0.0, 0.0), SVector(10.0, 0.0, 0.0)]
    @test ComplexMixtures.wrap!(x, x[2], box) ≈ SVector{3, Float64}[[11.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
    # Triclinic box
    box = CellListMap.Box(@SMatrix[10.0 5.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0], 1.0)
    x = [ SVector(1.0, 0.0, 0.0), SVector(10.0, 0.0, 0.0)]
    @test ComplexMixtures.wrap!(x, x[1], box) ≈ SVector{3, Float64}[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]] 
    x = [ SVector(1.0, 0.0, 0.0), SVector(10.0, 0.0, 0.0)]
    @test ComplexMixtures.wrap!(x, x[2], box) ≈ SVector{3, Float64}[[11.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
end

"""
  move!(x::AbstractVector, newcm::AbstractVector,beta, gamma, theta)

$(INTERNAL)

Translates and rotates a molecule according to the desired input center of coordinates and Euler rotations modifyies the vector x.

"""
function move!(x::AbstractVector{T}, newcm::AbstractVector{T}, beta, gamma, theta) where {T<:SVector}
    cm = mean(x)
    A = eulermat(beta, gamma, theta)
    for i in eachindex(x)
        x[i] = A * (x[i] - cm) + newcm
    end
    return nothing
end

@testitem "move!" begin
    using ComplexMixtures
    using StaticArrays

end

"""
  random_move!(x_ref::AbstractVector{T}, 
               irefatom::Int,
               sides::T,
               x_new::AbstractVector{T}, RNG) where {T<:SVector}

$(INTERNAL)

Function that generates a new random position for a molecule.

The new position is returned in `x_new`, a previously allocated array.

"""
function random_move!(x_ref::AbstractVector{T}, irefatom::Int, box::CellListMap.Box, x_new::AbstractVector{T}, RNG) where {T<:SVector}
    # To avoid boundary problems, the center of coordinates are generated in a 
    # much larger region, and wrapped aftwerwards
    scale = 100.0

    # Generate random coordiantes for the center of mass
    newcm = T(scale * (-sides[i] / 2 + random(RNG, Float64) * sides[i]) for i in 1:3)

    # Generate random rotation angles 
    beta = 2π * random(RNG, Float64)
    gamma = 2π * random(RNG, Float64)
    theta = 2π * random(RNG, Float64)

    # Copy the coordinates of the molecule chosen to the random-coordinates vector
    @. x_new = x_ref

    # Take care that this molecule is not split by periodic boundary conditions, by
    # wrapping its coordinates around its reference atom
    wrap!(x_new, x_ref[irefatom], box)

    # Move molecule to new position
    move!(x_new, newcm, beta, gamma, theta)

    return nothing
end



