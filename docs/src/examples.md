# [Examples](@id examples)

## List of examples

- [Protein in water/glycerol](@ref)
- [Polyacrylamide in DMDF](@ref)
- [POPC membrane in water/ethanol](@ref)
- [Glycerol/water mixture](@ref)

## How to run these examples

- Download and install [Julia](https://julialang.org)

- Install the required packages. Within Julia, do:
```julia-repl
julia> import Pkg

julia> Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings", "Formatting"])
```

- Get the files:
```bash
git clone https://github.com/m3g/ComplexMixturesExamples
```
The files associated to the following examples are distributed at [this page](https://github.com/m3g/ComplexMixturesExamples/tree/main/Protein_in_Glycerol). 

## Data

The [Data](https://github.com/m3g/ComplexMixturesExamples/tree/main/Protein_in_Glycerol/Data) directory contains the a pdb file of the system (`system.pdb`) and a sample from the trajectory (`glyc50.dcd`), with a few frames. It also contains the result of running the `mddf` calculation on the complete trajectory, `results_glyc50.json`. This last file was produced by `ComplexMixtures`, as indicated in the following examples. 

The sample trajectory is provided so that the first example can be run, yet do not expect that the results are the same, as the sampling is much lower in this case. The complete trajectory can be retrieved from [this link](https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing) (3GB file). 
