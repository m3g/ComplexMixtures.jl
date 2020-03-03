# Function that corrects the coordinates of the atoms in phantom
# cells for computation

function movephantomcoor(data :: DistanceData, iatom, ibox, jbox, kbox)

  if ibox == 0
    x = data.frame.x[iatom] - data.frame.sides[1]
  elseif ibox == data.lists.nboxes[1]+1
    x = data.frame.x[iatom] + data.frame.sides[1] 
  else 
    x = data.frame.x[iatom]
  end

  if jbox == 0
    y = data.frame.y[iatom] - data.frame.sides[2]
  elseif jbox == data.lists.nboxes[2]+1
    y = data.frame.y[iatom] + data.frame.sides[2]
  else
    y = data.frame.y[iatom]
  end

  if kbox == 0
    z = data.frame.z[iatom] - data.frame.sides[3]
  elseif kbox == data.lists.nboxes[3]+1
    z = data.frame.z[iatom] + data.frame.sides[3]
  else
    z = data.frame.z[iatom]
  end if

  return x, y, z

end

