#
# function that outputs the volumes and densities in the most
# natural units
#

@with_kw_noshow mutable struct Overview
  R :: Result
  domain_molar_volume :: Float64 = 0.
  density :: Density = Density()
  solvent_molar_volume :: Float64 = 0.
  solvent_molar_volume_bulk :: Float64 = 0.
  solute_molar_volume :: Float64 = 0.
end

function Base.show( io :: IO, ov :: Overview )
  println(bars)
  println(" Overview: ")
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
  println(" Simulation molar volume: $(ov.solute_molar_volume) cm³ mol⁻¹")
  println()
  println(" Using with dbulk=$(ov.R.dbulk) Å: ") 
  println(" Molar volume of the solute domain: $(ov.domain_molar_volume) cm³ mol⁻¹") 
  println()
  println(bars)
end

function overview(R :: Result)

  ov = Overview(R = R)

  ov.domain_molar_volume = R.volume.domain * units.Angs3tocm3permol

  ov.density.solute = R.density.solute * units.SitesperAngs3tomolperL 
  ov.density.solvent = R.density.solvent * units.SitesperAngs3tomolperL  
  ov.density.solvent_bulk = R.density.solvent_bulk * units.SitesperAngs3tomolperL   

  ov.solvent_molar_volume = 1000 / ov.density.solvent
  ov.solvent_molar_volume_bulk = 1000 / ov.density.solvent_bulk

  ov.solute_molar_volume = 1000 / ov.density.solute

  return ov
end
