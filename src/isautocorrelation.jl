#     
# just check if solute and solvent are the same
#
function isautocorrelation(trajectory :: Trajectory)
  if trajectory.solute.index == trajectory.solvent.index
    true
  else
    false
  end
end
