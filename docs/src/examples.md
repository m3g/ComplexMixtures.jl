# [Examples](@id examples)

## List of examples

- [Protein in water/glycerol](@ref)
- [Polyacrylamide in DMDF](@ref)
- [POPC membrane in water/ethanol](@ref)
- [Glycerol/water mixture](@ref)

## How to run these examples

1 Download and install [Julia](https://julialang.org)

To run the scripts, we suggest the following procedure:

2. Create a directory, for example `example1`.
3. Copy the required data files, indicated in each example.
4. Launch `julia` in that directory, activate the directory environment, and install the required packages. 
   This is done by launching Julia and executing:
   ```julia
   import Pkg 
   Pkg.activate(".")
   Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings, EasyFit"])
   exit()
   ```
5. Copy the code of each script in to a file, and execute with:
   ```julia
   julia -t auto script.jl
   ```
   Alternativelly (and perhaps preferrably), copy line by line the content of the script into
   the Julia REPL, to follow each step of the calculation. For a more advanced Julia usage,
   we suggest the [VSCode IDE](https://code.visualstudio.com/) with the 
   [Julia Language Support](https://www.julia-vscode.org/docs/dev/gettingstarted/) extension. 

