# [Updating scripts from v1 to v2](@id updating-scripts)

The syntax chances necessary to update script from version `1.X` to `2.X` of 
the package are:

## Atom selections

The previous `Selection` structure was renamed to `AtomSelection` for clarity.
- Before:
```julia
water = Selection(water; natomspermol=3)
```
- Now:
```julia
water = AtomSelection(water; natomspermol=3)
```

## Group contributions syntax

The syntax to computing group contributions is improved. Previusly, the `contrib` or
`contributions` functions required three somewhat redundant parameters. 
- Before:
The call to `contributions` required 3 parameters: the `Selection` structure,
the matrix of contributions, and the indexes of the atoms for which the
contributions were desired:
```julia
h_contributions = contributions(solvent, R.solvent_atom, h_indexes)
```
- Now:
The contributions are extracted from the `Result` data structure, by 
providing either a `SoluteGroup` or `SolventGroup` object, which are
setup with the group names, group indexes, atom names, or atom indexes:
```julia
h_contributions = contributions(R, SolventGroup(h_indexes))
```

## Frame weights
`frame_weights` is now an option of the `mddf` execution. That is previously,
they were defined in the `Options` data structure, and now they are passed
to the `mddf` function.
- Before:
```julia
options = Options(frame_weights=[1.0, 2.0])
results = mddf(trajectory, options)
```
- Now:
```julia
results = mddf(trajectory, options; frame_weights=[1.0, 2.0])
```