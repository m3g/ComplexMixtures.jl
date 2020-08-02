# [Results](@id results)

The results of a MDDF calculation are returned in a data structure which contains 
the MDDF, KB integrals, and atomic contributions. The following section
will assume that the computation was performed by calling the `mddf`
function with 
```julia
results = MDDF.mddf(trajectory)
``` 
such that the `results` variable contain the `Result` data structure. By
default, the histograms contain 500 bins (`binstep=0.002` and
`cutoff=10.`) such that all data-vectors will contain 500 lines.

## The Result data structure: main data

The most important data to be read from `resutls` are the distances,
minimum-distance distribution function, and KB integrals. These data is
stored in the following vectors: 

### Distances of the histograms: `results.d`
The following vector will contain values ranging from 0. to `cutoff`,
and the distance at each bin is the distance in that bin for which half
of the volume of the bin is within `d`, and half of the volume is above
`d`, if the volume was spherical: 
```julia
julia> results.d
500-element Array{Float64,1}:
 0.015874010519682
 0.033019272488946275
 ⋮
 9.970010030080179
 9.99001000999998

```

### Minimum-distance distribution function: `results.mddf`

The `results.mddf` vector will contain the main result, which the
minimum-distance distribution function. For a properly-sampled
simulation, it will be zero at very short distances and converge to 1.0
for distances smaller than the `cutoff`:
```julia
julia> results.mddf
500-element Array{Float64,1}:
 0.0
 0.0
     ⋮
 0.999052514965403
 1.001030818286187

```

A typical plot of `results.mddf` as a function of `results.d` will look
like:

```@raw html
<img src="../figures/mddf.png" width="60%">
```

Thus, this plot was obtained with the following code:
```julia
using Plots
plot(results.d,results.mddf,xlabel="d/A",ylabel="mddf(d) / L/mol") 
```

### Kirkwood-Buff integral: `results.kb`

The `results.kb` vector will contain the Kirkwood-Buff integral computed
as a function of the minimum-distance to the solute. For properly
sampled simulations, it is expected to converge at large distances.  
```julia
julia> results.kb
500-element Array{Float64,1}:
     0.0
    -4.3249356504752985
   -12.9804719721525
     ⋮
 35907.72186381783
 35860.13624162115

```

A typical plot of `results.kb` as a function of `results.d` will look
like:
```@raw html
<img src="../figures/kb.png" width="60%">
```

Thus, this plot was obtained with the following code:
```julia
using Plots
plot(results.d,results.kb,xlabel="d/A",ylabel="mddf(d) / L/mol") 
```

## Save and load results

Three functions serve the purpose of saving and loading the results
obtained with MDDF:

### Save data to recover it later 

```julia
MDDF.save(results,"results.json")
```
where `results` is the output data structure of the `MDDF.mddf()`
calculation, and `results.json` is the output file to be created. The
file is written in `JSON` format, thus is not naturally human-readable.

### Recover saved data

```julia
results = MDDF.read("results.json")
```
The `MDDF.read` function reads the output of the `save` function above,
and restores the results data structure.

### Write data to human-readable format

If you Want the results to be written as simple ASCII tables such that
you can read them with another analysis program, plotting graphic, or
just want to inspect the data visually, use:

```julia
MDDF.write(results,"results.dat")
```
Three files will be created by this function:

`results.dat`: Contains the main results, as the MDDF and KB-integral data.

`results-ATOM_CONTRIB_SOLVENT.dat`: contains the contribution of each
atom type of the solvent to the MDDF.

`results-ATOM_CONTRIB_SOLUTE.dat`: contains the contribution of each
atom type of the solute to the MDDF.







