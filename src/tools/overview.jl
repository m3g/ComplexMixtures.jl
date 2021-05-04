"""

Function that outputs the volumes and densities in the most natural units.


"""
@with_kw_noshow mutable struct Overview
  R::Result
  domain_molar_volume::Float64 = 0.
  density::Density = Density()
  solvent_molar_volume::Float64 = 0.
  solvent_molar_volume_bulk::Float64 = 0.
  solute_molar_volume::Float64 = 0.
end

function Base.show(io::IO, ov::Overview)
  println()
  println(bars)
  println()
  println(" MDDF Overview: ")
  println()
  println(" Solvent properties: ")
  println(" ------------------- ")
  println()
  println(" Simulation concentration: $(ov.density.solvent) mol L⁻¹")
  println(" Molar volume: $(ov.solvent_molar_volume) cm³ mol⁻¹")
  println()
  println(" Concentration in bulk: $(ov.density.solvent_bulk) mol L⁻¹")
  println(" Molar volume in bulk: $(ov.solvent_molar_volume_bulk) cm³ mol⁻¹ ")
  println()
  println(" Solute properties: ")
  println(" ------------------ ")
  println()
  println(" Simulation Concentration: $(ov.density.solute) mol L⁻¹")
  println(" Estimated solute partial molar volume: $(ov.solute_molar_volume) cm³ mol⁻¹")
  println()
  println(" Using with dbulk = $(ov.R.dbulk)Å: ") 
  println(" Molar volume of the solute domain: $(ov.domain_molar_volume) cm³ mol⁻¹") 
  println()
  println(" Auto-correlation: ",ov.R.autocorrelation)
  println()
  println(" Trajectory files and weights: ")
  for i in 1:length(ov.R.files)
    println("   $(ov.R.files[i]) - w = $(ov.R.weights[i])")
  end
  println()
  ifar = trunc(Int,ov.R.nbins - 1.0/ov.R.options.binstep)
  long_range_mean = mean(ov.R.mddf[ifar:ov.R.nbins])
  long_range_std = std(ov.R.mddf[ifar:ov.R.nbins])
  println(" Long range MDDF mean (expected 1.0): ", long_range_mean, " +/- ", long_range_std)
  long_range_mean = mean(ov.R.rdf[ifar:ov.R.nbins])
  long_range_std = std(ov.R.rdf[ifar:ov.R.nbins])
  println(" Long range RDF mean (expected 1.0): ", long_range_mean, " +/- ", long_range_std)
  println()
  println(bars)
end

function overview(R::Result)

  ov = Overview(R = R)

  # Molar volume of the solute domain
  ov.domain_molar_volume = R.volume.domain * units.Angs3tocm3permol

  # Density of the solute and of the solvent 
  ov.density.solute = R.density.solute * units.SitesperAngs3tomolperL 
  ov.density.solvent = R.density.solvent * units.SitesperAngs3tomolperL  
  ov.density.solvent_bulk = R.density.solvent_bulk * units.SitesperAngs3tomolperL   

  # Solvent molar volume
  ov.solvent_molar_volume = 1000 / ov.density.solvent
  ov.solvent_molar_volume_bulk = 1000 / ov.density.solvent_bulk

  # Solute molar volume computed from solvent density in bulk
  if R.autocorrelation
    ov.solute_molar_volume = ov.solvent_molar_volume
  else
    ov.solute_molar_volume = units.Angs3tocm3permol*
       (R.density.solvent_bulk*R.volume.total - R.solvent.nmols)/R.density.solvent_bulk
  end

  return ov
end



