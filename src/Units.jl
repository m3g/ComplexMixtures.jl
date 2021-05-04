"""

Unit conversions.

"""
@with_kw struct Units
  mole = 6.022140857e23 
  Angs3tocm3 = 1e24
  Angs3toL = 1e27
  Angs3tocm3permol = mole/Angs3tocm3
  Angs3toLpermol = mole/Angs3toL
  SitesperAngs3tomolperL = Angs3toL/mole
end
units = Units()

