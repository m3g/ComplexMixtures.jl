# [Results](@id results)

The results of a MDDF calculation are returned in a data structure which contains 
the MDDF, KB integrals, and atomic contributions. The following section
will assume that the computation was performed by calling the `mddf`
function with 
```julia
results = ComplexMixtures.mddf(trajectory)
``` 
such that the `results` variable contain the `Result` data structure. By
default, the histograms contain 500 bins (`binstep=0.002` and
`cutoff=10.`) such that all data-vectors will contain 500 lines.

To know how to save and load saved data, read the [next](@ref save) section.

## The Result data structure: main data

The most important data to be read from `resutls` are the distances,
minimum-distance distribution function, and KB integrals. These data is
stored in the following vectors:

## Distances of the histograms: `results.d`
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

## Minimum-distance distribution function: `results.mddf`

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

## Kirkwood-Buff integral: `results.kb`

The `results.kb` vector will contain the Kirkwood-Buff integral computed
as a function of the minimum-distance to the solute. For properly
sampled simulations, it is expected to converge at large distances.  
```julia
julia> results.kb
500-element Array{Float64,1}:
     0.0
    -0.3249356504752985
    -2.9804719721525
     ⋮
    0.72186381783
    1.13624162115

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

## Units

* The distance is assumed to be in ``\textrm{\AA}``, as this is the most common
  distance units in molecular simulations. The coordinates of the atoms
  are assumed be provided in ``\textrm{\AA}``. 

* The minimum-distance distribution function is unit-less, since it is the
  ratio of the density at each distance divided by an ideal-gas density.

* The Kirkwood-Buff integrals are returned in ``\textrm{cm}^3~\textrm{mol}^{-1}``, if the
  coordinates were provided in ``\textrm{\AA}``.

!!! warning
    If the coordinates are not in ``\textrm{\AA}``, the calculation will 
    proceed normaly, but the units of the KB integrals, which has units
    of volume per mol, should be
    converted to conform the length unit provided. 




