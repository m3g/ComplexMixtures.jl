
# Quick Guide

Of course, follow the [installation](@ref Installation) instructions first. 
Start `julia` and load the MDDF package, using:

```julia
using MDDF
```
And here we will use the `PDBTools` package to obtain the selections of
the solute and solvent molecules: 
```julia
using PDBTools
```
(see [Set solute and solvent](@ref selections) for details).

The fastest way to understand how to use this package is through an
example.  

Let us consider a system consisting of three components: a protein,
water, a cossolvent, and TMAO (trimetylamine-N-oxyde), which is a common
osmolyte known to stabilize protein structures.  

We want to study the interactions of the protein with TMAO, for example.
The computation of the MDDF is performed by defining the solute and
solvent selections, and running the calculation on the trajectory:

### Extract the indexes of the protein and TMAO atoms

The computation of the MDDF requires the indexes of the atoms of the
simulation which belong to the solute and the solvent, in this case the
protein and TMAO molecules. Here we will use the `PDBTools` package,
which provides a simple selection syntax.
```julia
protein_indexes = PDBTools.selindex("system.pdb","protein")
TMAO_indexes = PDBTools.selindex("system.pdb","resname TMAO")
```
The `protein_indexes` and `TMAO_indexes` are simply vectors containing
the indexes of the atoms of each selection (in a one-based scheme, that
is, the first atom is atom 1). 

!!! warning
    
    All the indexes in MDDF are 1-based. That means that the first
    atom in your structure file has index 1 in the coordinates vector.
    Please be careful when defining the selections.

### Define the protein as the solute

To define the protein solute, we need to provide a list of the indexes of the
atoms of the protein, and the number of protein molecules in the
simulation (most commonly 1):

```julia
solute = MDDF.Selection(protein_indexes,nmols=1)
```

### Define TMAO the solvent to be considered

Equivalently, the solvent is set up with:
```julia
solvent = MDDF.Selection(TMAO_indexes,natomspermol=14)
```
We need to provide either the number of TMAO molecules (with the `nmols`
keyword) or the number of atoms per molecule, as shown above. Since
there are typically many solvent molecules, it is easier to provide the
number of atoms per molecule with `natomspermol=14`. 

### Set the Trajectory structure

The `solute` and `solvent` data structures are then fed into the
`Trajectory` data structure, together with the trajectory file name,
with:
```julia
trajectory = MDDF.Trajectory("trajectory.dcd",solute,solvent)
```
In the case, the trajectory is of NAMD "dcd" format. All formats
supported by `Chemfiles` are automatically recognized. 

### Finally, run the computation and get the results:

If default options are used (as the bin size of the histograms, read all
frames without skipping any), just run the `mddf` with:
```julia
results = MDDF.mddf(trajectory)
```

### The results data structure obtained

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
distance to the protein.

The Kirkwood-Buff integral corresponding to that distribution is
provided in the `results.kb` vector.  

### Save the results

The results can be saved into a file (with JSON format) with:
```julia
MDDF.save(results,"./results.json")
```
And these results can be read aftwerwards with:
```julia
MDDF.read("./results.json")
```
Alternatively, a human-readable set of output files can be obtained to
be analyzed in other software (or plotted with alternative tools), with
```julia
MDDF.write(results,"./results.dat")
```

### Summary

The complete running example, therefore, should be:

```julia
# Load packages
using PDBTools
using MDDF 
using Plots

# Select indexes
protein_indexes = PDBTools.selindex("system.pdb","protein")
TMAO_indexes = PDBTools.selindex("system.pdb","resname TMAO")

# Setup solute and solvent structures
solute = MDDF.Selection(protein_indexes,nmols=1)
solvent = MDDF.Selection(TMAO_indexes,natomspermol=14)

# Setup the Trajectory structure
trajectory = MDDF.Trajectory("trajectory.dcd",solute,solvent)

# Run the calculation and get results
results = MDDF.mddf(trajectory)

# Save the reults to recover them later if required
MDDF.save(results,"./results.json")

# Plot the some of the most important results 
plot(results.d,results.mddf,xlabel="d",ylabel="MDDF") # plot the MDDF
savefig("./mddf.pdf")
plot(results.d,results.kb,xlabel="d",ylabel="KB") # plot the KB 
savefig("./kb.pdf")
```









