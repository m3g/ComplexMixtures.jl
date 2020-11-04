#
# Structure to contain trajectories read by the Chemfiles package. Many formats should
# be available through this interface, including the NamdDCD which are provided independently
# as an example. 
#

import Chemfiles

struct ChemFile{T<:Vf3} <: Trajectory

  #
  # Mandatory data for things to work
  #
  filename :: String
  format :: AbstractString
  stream :: Vector{Chemfiles.Trajectory} # mutable such that we can close it and open it again
  nframes :: Int64 

  # Vector wthat will be fed with the sizes of the periodic box 
  sides :: Vector{T}

  # Solute and solvent data
  solute :: Selection
  solvent :: Selection

  # Coordinates of the solute and solvent atoms in a frame (natoms,3) for each array:
  x_solute :: Vector{T}  # solute.natoms vectors of length 3
  x_solvent :: Vector{T} # solvent.natoms vectors of lenght 3

  #
  # Additional data required for input/output functions
  #
  natoms :: Int64

end

#
# Function open will set up the IO stream of the trajectory, fill up the 
# number of frames field and additional parameters if required 
#
# The IO stream must be returned in a position such that the "nextframe!" function
# will be able to read the first frame of the trajectory
#

function ChemFile( filename :: String, solute :: Selection, solvent :: Selection; format="" )

  stream = Vector{Chemfiles.Trajectory}(undef,1)
  stream[1] = Chemfiles.Trajectory(filename,'r',format)
  
  # Get the number of frames (the output of Chemfiles comes in UInt64 format, which is converted
  # to Int using (UInt % Int)
  nframes = Chemfiles.length(stream[1]) % Int
  
  # Read the first frame to get the number of atoms
  frame = Chemfiles.read(stream[1])
  natoms = Chemfiles.size(frame) % Int
  sides = Chemfiles.lengths(Chemfiles.UnitCell(frame)) # read the sides of the first frame
  Chemfiles.close(stream[1])

  # Reopen the stream, so that nextrame can read the first frame
  stream[1] = Chemfiles.Trajectory(filename,'r',format)

  return ChemFile( filename, # trajectory file name 
                   format, # trajectory format, is provided by the user
                   stream,
                   nframes, 
                   [ Vf3(sides) ], # array containing box sides
                   solute, solvent,
                   zeros(Vf3,solute.natoms),    
                   zeros(Vf3,solvent.natoms),  
                   natoms, # Total number of atoms
                 )
end

function Base.show( io :: IO, trajectory :: ChemFile )
  println(" ")
  println(" Trajectory read by Chemfiles with: ")
  println("    $(trajectory.nframes) frames.")
  println("    $(trajectory.natoms) atoms.")
  println("    PBC sides in current frame: $(trajectory.sides[1])")
end

#
# Function that reads the coordinates of the solute and solvent atoms from
# the next frame of the trajectory file 
#
# The function modifies sides, x_solute and x_solvent within the trajectory structure.
# Having these vectors inside the trajectory structure avoids having to allocate
# them everytime a new frame is read
#

function nextframe!( trajectory :: ChemFile ) 

  frame = Chemfiles.read(trajectory.stream[1])
  positions = Chemfiles.positions(frame)
  sides = Chemfiles.lengths(Chemfiles.UnitCell(frame))
  trajectory.sides[1] = Vf3(sides)

  # Save coordinates of solute and solvent in trajectory arrays (of course this could be avoided,
  # but the code in general is more clear aftwerwards by doing this)
  for i in 1:trajectory.solute.natoms
    trajectory.x_solute[i] = Vf3(positions[1,trajectory.solute.index[i]],
                                 positions[2,trajectory.solute.index[i]],
                                 positions[3,trajectory.solute.index[i]])
  end
  for i in 1:trajectory.solvent.natoms
    trajectory.x_solvent[i] = Vf3(positions[1,trajectory.solvent.index[i]],
                                  positions[2,trajectory.solvent.index[i]],
                                  positions[3,trajectory.solvent.index[i]])
  end

end

# Function that returns the sides of the periodic box given the data structure. In this
# case, returns the 3-element vector corresponding to the box sides of the given frame,
# it expected that the "nextframe" function fed this information already to the
# trajectory.sides vector in the current frame
getsides(trajectory :: ChemFile, iframe) = trajectory.sides[1]

#
# Function that closes the IO Stream of the trajectory
#
closetraj( trajectory :: ChemFile ) = Chemfiles.close( trajectory.stream[1] )

#
# Function that returns the trajectory in position to read the first frame
#
function firstframe( trajectory :: ChemFile )  
  Chemfiles.close(trajectory.stream[1])
  trajectory.stream[1] = Chemfiles.Trajectory(trajectory.filename,'r',trajectory.format)
end

