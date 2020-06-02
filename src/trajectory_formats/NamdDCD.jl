#
# Structure to contain DCD trajectories produces with Namd. 
#
# Must be mutable such that nframes, vector sizes, and additional parameters can be updated
# when the trajectory is first open
#

using FortranFiles

struct NamdDCD

  #
  # Mandatory data for things to work
  #
  filename :: String
  stream :: FortranFile # special type of stream required for reading DCD files
  nframes :: Int64 

  # This vector must be filled up with the size of the periodic cell, if it
  # is not defined in the DCD file. 
  sides :: Array{Float64}

  # Data structures of the solute and solvent 
  solute :: SoluteOrSolvent
  solvent :: SoluteOrSolvent

  # Coordinates of the solute and solvent atoms in a frame (natoms,3) for each array:
  x_solute :: Array{Float64}
  x_solvent :: Array{Float64}

  #
  # Additional properties that might be required for implementing IO (not necessary for every
  # input/output format)
  #
  sides_in_dcd :: Bool # if the DCD contains, or not, periodic cell information for each frame
  lastatom :: Int64 # The last atom to be read from each line

end

# This function initializes the structure above, returning the data and the vectors with
# appropriate lengths and, importantly, with the i/o stream OPENED, ready to read the first
# trajectory frame using the "nextframe" function.

function NamdDCD( filename :: String, solute :: SoluteOrSolvent, solvent :: SoluteOrSolvent )

  FortranDCD = FortranFile(filename)

  # Read header
  IntVec = Vector{Int32}(undef,17)
  hdr, read_nframes, IntVec[1:8], ts, IntVec[9:17] = read(FortranDCD, FString{4}, Int32, (Int32,8), Float64, (Int32,9))
  dummyi, title = read(FortranDCD, Int32, FString{80})
  read_natoms = read(FortranDCD,Int32)

  # Check if dcd file contains axis information
  sides_in_dcd = false
  x = 0.
  try
    x = read(FortranDCD, [ Float32 for i in 1:read_natoms ])
  catch err
    sides_in_dcd = true
  end

  # rewind and let it ready to read first frame in the first call to nextframe
  firstframe(FortranDCD)

  nframes = getnframes(FortranDCD,sides_in_dcd) 
  stream = FortranDCD
  lastatom = max(maximum(solute.index),maximum(solvent.index))

  return NamdDCD( filename, stream, nframes, 
                  Vector{Float64}(undef,3), # sides
                  solute, solvent,
                  Array{Float64}(undef,solute.natoms,3), # solute atom coordinates
                  Array{Float64}(undef,solvent.natoms,3), # solvent atom coordinates
                  sides_in_dcd, lastatom )

end

import Base.show
function Base.show( io :: IO, traj :: NamdDCD )
  println(" Trajectory of NamdDCD format containing: ")
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

function nextframe!( trajectory:: NamdDCD )

  # Read the sides of the box from the DCD file, otherwise they must be set manually before
  if trajectory.sides_in_dcd
    sides_read = read(trajectory.stream,(Float64,6))
    trajectory.sides[1] = sides_read[1]
    trajectory.sides[2] = sides_read[3]
    trajectory.sides[3] = sides_read[6]
  end
  
  # Read the coordinates  
  x = read(trajectory.stream,(Float32,trajectory.lastatom))
  y = read(trajectory.stream,(Float32,trajectory.lastatom))
  z = read(trajectory.stream,(Float32,trajectory.lastatom))

  # Save coordinates of solute and solvent in trajectory arrays
  for i in 1:solute.natoms
    trajectory.x_solute[i,1] = x[trajectory.solute.index[i]]
    trajectory.x_solute[i,2] = y[trajectory.solute.index[i]]
    trajectory.x_solute[i,3] = z[trajectory.solute.index[i]]
  end
  for i in 1:solvent.natoms
    trajectory.x_solvent[i,1] = x[trajectory.solvent.index[i]]
    trajectory.x_solvent[i,2] = y[trajectory.solvent.index[i]]
    trajectory.x_solvent[i,3] = z[trajectory.solvent.index[i]]
  end

end

#
# Function that closes the IO Stream of the trajectory
#

function close( trajectory :: NamdDCD )
  close(trajectory.stream)
end

# Function that returns the sides of the periodic box given the data structure
# In this case, just return the sides vector which 

function getsides(trajectory :: NamdDCD, iframe)
  # In this (most common) case, sides is a vector and must only be returned
  if trajectory.sides_in_dcd
    return trajectory.sides
  # otherwise, sides is an array that contains the sides for each frame, and we return the
  # vector containing the sides of the current fraem
  else
    return [ trajectory.sides[iframe,1], trajectory.sides[iframe,2], trajectory.sides[iframe,3] ]
  end
end

#
# Auxiliary functions
#

#
# Sometimes the DCD files contains a wrong number of frames in the header, so to
# get the actual number of frames, it is better to read it
#

function getnframes(FortranDCD :: FortranFile, sides_in_dcd :: Bool )
  firstframe(FortranDCD)
  nframes = 0
  while true
    try 
      if sides_in_dcd
        read(FortranDCD,Float64)
      end
      read(FortranDCD,Float32)
      read(FortranDCD,Float32)
      read(FortranDCD,Float32)
      nframes = nframes + 1
    catch
      firstframe(FortranDCD)
      return nframes
    end
  end
end

#
# Leave DCD file in position to read the first frame
#
function firstframe(FortranDCD :: FortranFile)
    # rewind
    rewind(FortranDCD)
    # skip header
    read(FortranDCD)
    read(FortranDCD)
    read(FortranDCD)
end



