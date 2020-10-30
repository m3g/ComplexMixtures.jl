#
# Structure to contain PDB trajectories. Frames separated by "END", and with periodic cell sizes
# in the "CRYST1" field.
#
# These structure and functions can be used as a template to implement the reading of other 
# trajectory formats
#

struct PDBTraj <: Trajectory

  #
  # Mandatory data for things to work
  #
  filename :: String
  stream :: IOStream
  nframes :: Int64 

  # This Array (3,nframes) must be filled up with the size of the periodic cell
  sides :: Array{Float64}

  # Solute and solvent data
  solute :: Selection
  solvent :: Selection

  # Coordinates of the solute and solvent atoms in a frame (natoms,3) for each array:
  x_solute :: Array{Float64,2}  # (3,solute.natoms)
  x_solvent :: Array{Float64,2} # (3,solvent.natoms)

  #
  # Additional data required for input/output functions
  #
  natoms :: Int64

  # Auxiliary vectors to contain all coordinates of a frame on reading 
  x_read :: Array{Float64,2} 

end

#
# Function open will set up the IO stream of the trajectory, fill up the 
# number of frames field and additional parameters if required 
#
# The IO stream must be returned in a position such that the "nextframe!" function
# will be able to read the first frame of the trajectory
#

function PDBTraj( pdbfile :: String, solute :: Selection, solvent :: Selection)

  stream = open(pdbfile,"r")
  
  # Get the number of frames and number of atoms, check for adequacy of the format
  nframes = 0
  natoms = 0
  for line in eachline(stream)
    if line[1:3] == "END"
      nframes = nframes + 1
    elseif nframes == 0 && (line[1:4] == "ATOM" || line[1:6] == "HETATM")
      natoms = natoms + 1
    end
  end
  close(stream)

  # Fill-up the sides vector of the trajectory. We assume here
  # that the sides are stored for each frame in the "CRYST1" field, for each frame.
  # Here, we exemplify the option to read all sides at once and store them in a 
  # sides(nframes,3) array. Alternatively, the sides could be read for each frame
  # independently within the "nextframe!" function, and saved as a sides(3) vector.  
  # The function "getsides", below, must be adapted accordingly to return the correct
  # sides of the periodic box in each frame.

  sides = zeros(Float64,3,nframes)
  stream = open(pdbfile,"r")
  iframe = 0
  for line in eachline(stream)
    s = split(line)
    if s[1] == "CRYST1"
      iframe = iframe + 1
      sides[1,iframe] = parse(Float64,s[2])
      sides[2,iframe] = parse(Float64,s[3])
      sides[3,iframe] = parse(Float64,s[4])
    end
  end
  close(stream)

  # Open again the stream, so that nextrame can read the first frame
  stream = open(pdbfile,"r")

  return PDBTraj( pdbfile, # trajectory file name 
                  stream,
                  nframes, 
                  sides, # array containing box sides
                  solute, solvent,
                  zeros(Float64,3,solute.natoms),    
                  zeros(Float64,3,solvent.natoms),  
                  natoms, # Total number of atoms
                  Array{Float64}(undef,3,natoms) # Auxiliary array for reading
                )
end

function Base.show( io :: IO, trajectory :: PDBTraj )
  println(" Trajectory in PDB format with: ")
  println("    $(trajectory.nframes) frames.")
  println("    $(trajectory.natoms) atoms.")
end

#
# Function that reads the coordinates of the solute and solvent atoms from
# the next frame of the trajectory file 
#
# The function modifies sides, x_solute and x_solvent within the trajectory structure.
# Having these vectors inside the trajectory structure avoids having to allocate
# them everytime a new frame is read
#

function nextframe!( trajectory :: PDBTraj ) 

  iatom = 0
  line = [ "START" ]
  while line[1] != "END"
    record = readline(trajectory.stream)
    line = split(record)
    if line[1] == "ATOM" || line[1] == "HETATM"
      iatom = iatom + 1
      trajectory.x_read[1,iatom] = parse(Float64,record[31:38])
      trajectory.x_read[2,iatom] = parse(Float64,record[39:46])
      trajectory.x_read[3,iatom] = parse(Float64,record[47:54])
    end
  end

  # Save coordinates of solute and solvent in trajectory arrays
  for i in 1:trajectory.solute.natoms
    trajectory.x_solute[1,i] = trajectory.x_read[1,trajectory.solute.index[i]]
    trajectory.x_solute[2,i] = trajectory.x_read[2,trajectory.solute.index[i]]
    trajectory.x_solute[3,i] = trajectory.x_read[3,trajectory.solute.index[i]]
  end
  for i in 1:trajectory.solvent.natoms
    trajectory.x_solvent[1,i] = trajectory.x_read[1,trajectory.solvent.index[i]]
    trajectory.x_solvent[2,i] = trajectory.x_read[2,trajectory.solvent.index[i]]
    trajectory.x_solvent[3,i] = trajectory.x_read[3,trajectory.solvent.index[i]]
  end

end

# Function that returns the sides of the periodic box given the data structure. In this
# case, returns the 3-element vector corresponding to the box sides of the given frame

function getsides(trajectory :: PDBTraj, iframe)
  # Sides is expected to be an array that contains the sides for each frame, and we return the
  # vector containing the sides of the current fraem
  return @view(trajectory.sides[1:3,iframe]) 
end

#
# Function that closes the IO Stream of the trajectory
#
closetraj( trajectory :: PDBTraj ) = close( trajectory.stream )

#
# Function that returns the trajectory in position to read the first frame
#
firstframe( trajectory :: PDBTraj ) = seekstart( trajectory.stream )






