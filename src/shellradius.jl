# Compute the point in which the radius comprises half of the
# volume of the shell

function shellradius(i,step)

  rmin = (i-1)*step
  return ( 0.5*( (rmin+step)^3 + rmin^3 ) )**(1./3.)

end 

