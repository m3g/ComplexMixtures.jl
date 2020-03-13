#
# Structures to contain the results of the MDDF calculation
#

mutable struct Density

  solute :: Float64
  solvent :: Float64
  solvent_bulk :: Float64

end

mutable struct Volume

  shell :: Vector{Float64}
  total :: Float64
  bulk :: Float64

end

struct MDDF_Results

  d :: Vector{Float64}
  mddf :: Vector{Float64}
  kb :: Vector{Float64}

  solute_atom :: Array{Float64}
  solvent_atom :: Array{Float64}
  
  count :: Vector{Float64}
  count_random :: Vector{Float64}
  
  density :: Density
  volume :: Volume

end

# Initializer from nbins and the number of atoms of the solute and solvent molecules

function MDDF_Results( nbins :: Int64, solute :: Solute, solvent :: Solvent ) 
  return MDDF_results( zeros(Float64,nbins), # d
                       Vector{Float64}(undef,nbins), # mddf
                       Vector{Float64}(undef,nbins), # kb
                       zeros(Float64,solute.natomspermol,nbins), # solute_atom
                       zeros(Float64,solvent.natomspermol,nbins), # solvent_atom
                       zeros(Float64,nbins), # count
                       zeros(Float64,nbins), # count_random
                       Density( 0. , # solute
                                0. , # solvent
                                0.   # solent_bulk
                              ),
                       Volume( Vector{Float64}(undef,nbins),
                               0. , # total
                               0.   # bulk
                             )
                     )
end


