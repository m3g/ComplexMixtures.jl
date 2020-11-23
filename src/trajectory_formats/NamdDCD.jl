#
# Structure to contain DCD trajectories produces with Namd. 
#

import FortranFiles

struct NamdDCD{T<:AbstractVector} <: Trajectory

  #
  # Mandatory data for things to work
  #
  filename :: String
  stream :: FortranFile # special type of stream required for reading DCD files
  nframes :: Int64 

  # This vector must be filled up with the size of the periodic cell, if it
  # is not defined in the DCD file. 
  sides :: Vector{T}

  # Data structures of the solute and solvent 
  solute :: Selection
  solvent :: Selection

  # Coordinates of the solute and solvent atoms in a frame (3,natoms) for each array:
  x_solute :: Vector{T}
  x_solvent :: Vector{T}

  #
  # Additional properties that might be required for implementing IO (not necessary for every
  # input/output format)
  #
  sides_in_dcd :: Bool # if the DCD contains, or not, periodic cell information for each frame
  lastatom :: Int64 # The last atom to be read from each line

  # Auxiliary vectors to read the coordinates without having the allocate/deallocate every time
  sides_read :: Vector{Float64}
  x_read :: Vector{Float32}
  y_read :: Vector{Float32}
  z_read :: Vector{Float32}

end

# This function initializes the structure above, returning the data and the vectors with
# appropriate lengths and, importantly, with the i/o stream OPENED, ready to read the first
# trajectory frame using the "nextframe" function.

function NamdDCD( filename :: String, solute :: Selection, solvent :: Selection; 
                  T :: Type = SVector{3,Float64})

  stream = FortranFile(filename)

  # Read header
  IntVec = Vector{Int32}(undef,17)
  hdr, read_nframes, IntVec[1:8], ts, IntVec[9:17] = read(stream, FString{4}, Int32, (Int32,8), Float64, (Int32,9))
  dummyi, title = read(stream, Int32, FString{80})
  read_natoms = read(stream,Int32)

  # Check if dcd file contains axis information
  sides_in_dcd = false
  x = 0.
  try
    x = read(stream, [ Float32 for i in 1:read_natoms ])
  catch err
    sides_in_dcd = true
  end

  # rewind and let it ready to read first frame in the first call to nextframe
  firstframe(stream)

  nframes = getnframes(stream,sides_in_dcd) 
  lastatom = max(maximum(solute.index),maximum(solvent.index))

  # Most commonly the sides of the box are written in each frame of the DCD file, and will
  # be updated upon reading the frame. Alternatively, the user must provide the sides in all
  # frames by filling up an array with the box side data.
  if sides_in_dcd
    sides = zeros(T,1) 
  else
    sides = zeros(T,nframes) 
  end

  return NamdDCD( filename, stream, nframes, 
                  sides, # sides vector (if in dcd) or array to be filled up later
                  solute, solvent,
                  zeros(T,solute.natoms), # solute atom coordinates
                  zeros(T,solvent.natoms), # solvent atom coordinates
                  sides_in_dcd, lastatom,
                  Vector{Float64}(undef,6), # auxiliary vector to read sides
                  Vector{Float32}(undef,lastatom), # auxiliary x
                  Vector{Float32}(undef,lastatom), # auxiliary y
                  Vector{Float32}(undef,lastatom)  # auxiliary z
                )

end

function Base.show( io :: IO, traj :: NamdDCD )
  println(" Trajectory in NamdDCD format containing: ")
  println("     $(traj.nframes) frames ") 
  println("     Sides in DCD: $(traj.sides_in_dcd) ") 
end

#
# Function that reads the coordinates of the solute and solvent atoms from
# the next frame of the trajectory file 
#
# The function modifies sides, x_solute and x_solvent within the trajectory structure.
# Having these vectors inside the trajectory structure avoids having to allocate
# them everytime a new frame is read
#

function nextframe!(trajectory :: NamdDCD{T}) where T

  # Read the sides of the box from the DCD file, otherwise they must be set manually before
  if trajectory.sides_in_dcd
    read(trajectory.stream,trajectory.sides_read)
    trajectory.sides[1] = T(trajectory.sides_read[1],
                            trajectory.sides_read[3],
                            trajectory.sides_read[6])
  end
  
  # Read the coordinates  
  read(trajectory.stream,trajectory.x_read)
  read(trajectory.stream,trajectory.y_read)
  read(trajectory.stream,trajectory.z_read)

  # Save coordinates of solute and solvent in trajectory arrays
  for i in 1:trajectory.solute.natoms
    trajectory.x_solute[i] = T(trajectory.x_read[trajectory.solute.index[i]],
                               trajectory.y_read[trajectory.solute.index[i]],
                               trajectory.z_read[trajectory.solute.index[i]])
  end
  for i in 1:trajectory.solvent.natoms
    trajectory.x_solvent[i] = T(trajectory.x_read[trajectory.solvent.index[i]],
                                trajectory.y_read[trajectory.solvent.index[i]],
                                trajectory.z_read[trajectory.solvent.index[i]])
  end

  nothing
end

#
# Function that closes the IO Stream of the trajectory
#

function closetraj( trajectory :: NamdDCD )
  FortranFiles.close(trajectory.stream)
end

#
# Function that returns a vector of dimension 3 with the sides of the periodic box 
# given the way that the box side information is stored in the Trajectory structure
#

function getsides(trajectory :: NamdDCD, iframe)
  # In this (most common) case, sides is a vector and must only be returned
  if trajectory.sides_in_dcd
    return trajectory.sides[1]
  # otherwise, sides is an array that contains the sides for each frame, and we return the
  # vector containing the sides of the current frame
  else
    return trajectory.sides[iframe]
  end
end

#
# Leave DCD file in position to read the first frame
#

function firstframe(stream :: FortranFile)
    # rewind
    rewind(stream)
    # skip header
    read(stream)
    read(stream)
    read(stream)
end
firstframe( trajectory :: NamdDCD ) = firstframe( trajectory.stream )

#
# Auxiliary functions
#

#
# Sometimes the DCD files contains a wrong number of frames in the header, so to
# get the actual number of frames, it is better to read it
#

function getnframes(stream :: FortranFile, sides_in_dcd :: Bool )
  firstframe(stream)
  nframes = 0
  while true
    try 
      if sides_in_dcd
        read(stream,Float64)
      end
      read(stream,Float32)
      read(stream,Float32)
      read(stream,Float32)
      nframes = nframes + 1
    catch
      firstframe(stream)
      return nframes
    end
  end
end


