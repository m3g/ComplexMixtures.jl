#
# Structure to contain DCD trajectories produces with Namd. 
#
# Must be mutable such that nframes, vector sizes, and additional parameters can be updated
# when the trajectory is first open
#

mutable struct NamdDCD

  #
  # Mandatory data for things to work
  #
  filename :: String
  iostream :: IOStream
  nframes :: Int64 

  # This vector must be filled up with the size of the periodic cell, if it
  # is not defined in the DCD file. 
  sides :: Array{Float64}

  # Coordinates of the solute and solvent atoms in a frame (natoms,3) for each array:
  x_solute :: Array{Float64}
  x_solvent :: Array{Float64}

  #
  # Additional properties that might be required for implementing IO
  #
  sides_in_dcd :: Bool # if the DCD contains, or not, periodic cell information for each frame
  lastatom :: Int64 # The last atom to be read from each line

end

using FortranFiles

#
# Function open will set up the IO stream of the trajectory, fillup the 
# number of frames field and additional parameters if required (in this case,
# if the trajectory file has or not pbc information). 
#
# The IO stream must be returned in a position such that the "nextframe" function
# will be able to read the first frame of the trajectory
#

function open( trajectory :: NamdDCD )

  dcdfile = FortranFile(trajectory.filename)
  
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
  firstframe(dcdfile)

  # Modify the trajectory structure with the data read:
  trajectory.nframes = getnframes(dcdfile)
  trajectory.iostream = dcdfile
  trajectory.sides_in_dcd = sides_in_dcd
  
end

#
# Function that reads the coordinates of the solute and solvent atoms from
# the next frame of the trajectory file 
#
# The function modifies sides, x_solute and x_solvent within the trajectory structure.
# Having these vectors inside the trajectory structure avoids having to allocate
# them everytime a new frame is read
#

function nextframe!( trajectory:: NamdDCD, solute :: Solute, solvent :: Solvent )

  # Read the sides of the box from the DCD file, otherwise they must be set manually before
  if trajectory.sides_in_dcd
    sides_read = read(trajectory.iostream,(Float64,6))
    trajectory.sides = [ sides_read[1], sides_read[3], sides_read[6] ]
  end
  
  # Read the coordinates  
  x = read(dcdfile,(Float32,lastatom))
  y = read(dcdfile,(Float32,lastatom))
  z = read(dcdfile,(Float32,lastatom))

  # Save coordinates of solute and solvent in trajectory arrays
  for i in 1:solute.n
    trajectory.x_solute[i,1] = x[solute.index[i]]
    trajectory.x_solute[i,2] = y[solute.index[i]]
    trajectory.x_solute[i,3] = z[solute.index[i]]
  end
  for i in 1:solvent.n
    trajectory.x_solvent[i,1] = x[solvent.index[i]]
    trajectory.x_solvent[i,2] = y[solvent.index[i]]
    trajectory.x_solvent[i,3] = z[solvent.index[i]]
  end

end

#
# Function that closes the IO Stream of the trajectory
#

function close( trajectory :: NamdDCD )
  close(trajectory.iostream)
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

function getnframes(FortranDCD :: FortranFile, dcdaxis :: Bool )
  firstframe(FortranDCD)
  nframes = 0
  while true
    try 
      if dcdaxis
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



