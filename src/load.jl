"""

```
load(filename::String)
```

Function to load the json saved results file into the `Result` data structure.

"""
function load(filename::String)
  f = open(filename,"r")
  R = JSON3.read(f,MutableResult)
  # Need to reshape the solute and solvent atom contributions, because
  # the data is read in a single column
  solute_atom = reshape(R.solute_atom,R.nbins,:)
  solvent_atom = reshape(R.solvent_atom,R.nbins,:)
  R.solute_atom = solute_atom
  R.solvent_atom = solvent_atom
  return Result([getfield(R,field) for field in fieldnames(Result)]...)
end

