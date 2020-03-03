#
# Function that wraps the coordinates to obtain minimum images
# around a defined center
# 
# It modifies the x,y,z input vectors
#

function wrap!(sides,x,y,z;center=[0.,0.,0.],sel=[-1])

  # If no selection was set, wrap all
  if sel[1] == -1 
    n = length(x)
    for i in 1:n
      wrapone!(i,sides,x,y,z,center)
    end

  # Wrap only the selection
  else
    n = length(sel)
    for i in 1:n
      wrapone!(sel[i],sides,x,y,z,center)
    end
  end

end

function wrapone!(i,sides,x,y,z,center)

  x[i] = x[i] - center[1]
  y[i] = y[i] - center[2]
  z[i] = z[i] - center[3]

  x[i] = x[i]%sides[1]
  y[i] = y[i]%sides[2]
  z[i] = z[i]%sides[3]

  if x[i] > sides[1]/2 ; x[i] = x[i] - sides[1] ; end
  if y[i] > sides[2]/2 ; y[i] = y[i] - sides[2] ; end
  if z[i] > sides[3]/2 ; z[i] = z[i] - sides[3] ; end

  if x[i] < -sides[1]/2 ; x[i] = x[i] + sides[1] ; end
  if y[i] < -sides[2]/2 ; y[i] = y[i] + sides[2] ; end
  if z[i] < -sides[3]/2 ; z[i] = z[i] + sides[3] ; end

  x[i] = x[i] + center[1]
  y[i] = y[i] + center[2]
  z[i] = z[i] + center[3]

end
