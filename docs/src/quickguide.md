
# Quick Guide

Of course, follow the [installation](@ref Installation) instructions first. 
A complete working example is shown below, and in the section that follows each 
command is described in detail.

## Complete example

Here we show the input file required for the study of the solvation of a protein
by the `TMAO` solvent, which is a molecule 4 atoms. The protein is assumed to be
at infinite dilution in the simulation. The trajectory of the simulation is in `DCD`
format in this example, which is the default output of `NAMD` and `CHARMM` simulation
packages.

```julia
# Activate environment (see the Installation -> Recommended Workflow manual section)
import Pkg; Pkg.activate("/home/user/MyNewPaper")

# Load packages
using PDBTools
using ComplexMixtures 
using Plots

# Load PDB file of the system
atoms = readPDB("./system.pdb")

# Select the protein and the TMAO molecules
protein = select(atoms,"protein")
tmao = select(atoms,"resname TMAO")

# Setup solute and solvent structures
solute = Selection(protein,nmols=1)
solvent = Selection(tmao,natomspermol=14)

# Setup the Trajectory structure
trajectory = Trajectory("./trajectory.dcd",solute,solvent)

# Run the calculation and get results
results = mddf(trajectory)

# Save the reults to recover them later if required
save(results,"./results.json")

# Plot the some of the most important results 
plot(results.d,results.mddf,xlabel="d",ylabel="MDDF") # plot the MDDF
savefig("./mddf.pdf")
plot(results.d,results.kb,xlabel="d",ylabel="KB") # plot the KB 
savefig("./kb.pdf")
```

## Running the example

Given that this code is saved into a file named `example.jl`, 
it can be run within the Julia REPL with:
```julia
julia> include("example.jl")
```
or directly with:
```
julia -t auto example.jl
```
where `-t auto` will launch `julia` with multi-threading. It is highly recommended to use multi-threading!

!!! note
    Julia compiles the code the first time it is run in each section. Thus, running script from the command line with, for example, `julia -t auto example.jl` may appear slow, particularly if you are modifying the script interactively. Ideally, *do not restart Julia*, just repeatedly include your script in the same Julia section. The second time the script is loaded will be much faster. For example:
    ```julia-repl
    julia> @time include("./example.jl") 
    # ... some output
    27.554095 seconds (72.13 M allocations: 4.623 GiB, 4.37% gc time, 11.96% compilation time)
    julia> @time include("./example.jl")
    # ... some output
    0.703780 seconds (3.24 M allocations: 897.260 MiB, 12.22% gc time)
    ```
    The first time the script was called it took ~30s, which included compilation of the code and loading of the packages. The second time the script was included it took 0.7s. Thus, for interactive modification of the script, don't restart Julia.

    This has been much improved in Julia 1.9, particularly if a stable environment is used, as suggested.

## Detailed description of the example

Start `julia` and load the ComplexMixtures package, using:

```julia
using ComplexMixtures
```
And here we will use the `PDBTools` package to obtain the selections of
the solute and solvent molecules: 
```julia
using PDBTools
```
(see [Set solute and solvent](@ref selections) for details).

The fastest way to understand how to use this package is through an
example.  

Let us consider a system of three components: a protein,
water, a cossolvent: TMAO (trimetylamine-N-oxyde), which is a common
osmolyte known to stabilize protein structures. A picture of this system
is shown below, with the protein in blue, water, and TMAO molecules. The
system was constructed with [Packmol](http://m3g.iqm.unicamp.br/packmol)
and the figure was produced with
[VMD](https://www.ks.uiuc.edu/Research/vmd/).

```@raw html
<center>
<img src="../figures/proteinTMAO.png" width=60%>
</center>
```

We want to study the interactions of the protein with TMAO in this example.
The computation of the MDDF is performed by defining the solute and
solvent selections, and running the calculation on the trajectory.

### Define the protein as the solute

To define the protein as the solute, we will use the PDBTools package,
which provides a handy selection syntax. First, read the PDB file using 
```julia
atoms = readPDB("./system.pdb")
```
Then, let us select the protein atoms (here we are using the `PDBTools.select` function):
```julia
protein = select(atoms,"protein")
```
And, finally, let us use the `Selection` function to setup the
structure required by the MDDF calculation:
```julia
solute = Selection(protein,nmols=1)
```

!!! note
    It is necessary to indicate how many molecules (in this case,
    `nmols=1`, so that ComplexMixtures knows that the solute is to be considered
    as single structure. In this case there is no ambiguity, but if
    the solute was a micelle, for example, this option would let 
    ComplexMixtures know that one wants to consider the micelle as a single structure.


### Define TMAO the solvent to be considered

Equivalently, the solvent is set up with:
```julia
tmao = select(atoms,"resname TMAO")
solvent = Selection(tmao,natomspermol=14)

```

!!! note
    Here we opted to provide the number of atoms of a TMAO molecules (with the
    `natomspermol` keyword). This is generally more practical for small
    molecules than to provide the number of molecules.

### Set the Trajectory structure

The `solute` and `solvent` data structures are then fed into the
`Trajectory` data structure, together with the trajectory file name,
with:
```julia
trajectory = Trajectory("trajectory.dcd",solute,solvent)
```
In the case, the trajectory is of NAMD "dcd" format. All formats
supported by [Chemfiles](http://chemfiles.org/Chemfiles.jl/latest/) 
are automatically recognized. 

### Finally, run the computation and get the results:

If default options are used (as the bin size of the histograms, read all
frames without skipping any), just run the `mddf` with:
```julia
results = mddf(trajectory)

```
Some optional parameters for the computation are available in the
[Options](@ref options) section.

### The `results` data structure obtained

The `results` data structure contains all the results of the MDDF
calculation, including:

`results.d` : Vector containing the distances to the solute. 

`results.mddf` : Vector containing the minimum-distance distribution
function at each distance.

That means, for example, that 
```julia
plot(results.d,results.mddf,xlabel="d / \AA",ylabel="MDDF") 

```
results in the expected plot of the MDDF of TMAO as a function of the
distance to the protein:

```@raw html
<center>
<img src="../figures/mddf.png" width="60%">
</center>
```

The Kirkwood-Buff integral corresponding to that distribution is
provided in the `results.kb` vector, and can be also directly plotted 
with   

```julia
plot(results.d,results.kb,xlabel="d / \AA",ylabel="MDDF") 
```
to obtain:

```@raw html
<center>
<img src="../figures/kb.png" width="60%">
</center>
```

See the [Atomic and group contributions](@ref contrib) section for a
detailed account on how to obtain a molecular picture of the solvation
by splitting the MDDF in the contributions of each type of atom of the
solvent, each type of residue of the protein, etc.

### Save the results

The results can be saved into a file (with JSON format) with:
```julia
save(results,"./results.json")
```
And these results can be loaded afterwards with:
```julia
load("./results.json")
```
Alternatively, a human-readable set of output files can be obtained to
be analyzed in other software (or plotted with alternative tools), with
```julia
write(results,"./results.dat")
```








