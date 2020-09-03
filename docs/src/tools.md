# [Tools](@id Tools)

Some tools are provided to analyze the results:

## Computing radial distribution functions

The distributions returned by the `mddf` function (meaning `mddf` and
`rdf`) vectors, are normalized either the random reference state or
using a site count based on the numerical integration of the volume
corresponding to each minimum-distance to the solute. If, however, the
solute is defined by a single atom (as the oxygen atom of water, for
example), the numerical integration of the volume can be replaced by a
simple analytical spherical shell volume, reducing noise. The `gr`
function returns the radial distribution function and the KB integral 
computed from the results, using this volume estimate: 

```julia
g, kb = ComplexMixtures.gr(R)

```

By default, the single-reference count (`rdf_count`) of the Result
structure will be used to compute the radial distribution function. The
function can be called with explicit control of all input parameters: 

```julia
g, kb = ComplexMixtures.gr(r,count,density,binstep)

```
where:

| Parameter | Definition | Result structure output data to provide |
|:---------:|:-----------|:-------------------|
| `r`       | Vector of distances | The `d` vector |
| `count`   | Number of site counts at each r | The `rdf` or `mddf` vectors  |
| `density` | Bulk density | The `density.solvent_bulk` or `density.solvent` densities. |
| `binstep` | The histogram step | The `options.binstep`  |
|           |                    |                        |

Example:
```julia
...
R = ComplexMixtures.mddf(trajectory,options)
g, kb = ComplexMixtures.gr(R.d,R.rdf_count,R.density.solvent_bulk,R.options.binstep)

```




