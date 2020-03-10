# Computes the volume of the spherical shell 
# defined within [(i-1)*step,i*step]

function sphericalshellvolume(i,step)

  fourthirdsofpi = (4./3.)*pi
  rmin = (i-1)*step
  return fourthirdsofpi*( (rmin+step)^3 - rmin^3 )

end


