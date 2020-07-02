# Function to compute Euclidean distances

distance(x,y) = sqrt((y[1]-x[1])^2 + (y[2]-x[2])^2 + (y[3]-x[3])^2)

distance(x,y,j) = sqrt((y[j,1]-x[1])^2 + (y[j,2]-x[2])^2 + (y[j,3]-x[3])^2)

distance(x,y,i,j) = sqrt((y[j,1]-x[i,1])^2 + (y[j,2]-x[i,2])^2 + (y[j,3]-x[i,3])^2)

# With possible periodic conditions 

function distance( box :: Box, x, y )
 
  dx = (x[1]-y[1])%box.sides[1]
  dy = (x[2]-y[2])%box.sides[2]
  dz = (x[3]-y[3])%box.sides[3]

  (dx > box.sides[1]/2) && (dx = dx - box.sides[1])
  (dy > box.sides[2]/2) && (dx = dx - box.sides[2])
  (dz > box.sides[3]/2) && (dx = dx - box.sides[3])

  (dx < -box.sides[1]/2) && (dx = dx + box.sides[1])
  (dy < -box.sides[2]/2) && (dx = dx + box.sides[2])
  (dz < -box.sides[3]/2) && (dx = dx + box.sides[3])

  return sqrt(dx^2+dy^2+dz^2)

end




