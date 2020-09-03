#
# Unit conversion
#
@with_kw struct Units
  mole = 6.022140857e23 
  Angs3tocm3 = 1e24
  Angs3tocm3permol = mole/Angs3tocm3
end
units = Units()

