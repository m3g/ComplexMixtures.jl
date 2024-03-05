
# Quick Guide

Of course, follow the [installation](@ref Installation) instructions first. 
A complete working example is shown below, and in the section that follows each 
command is described in detail.

## Basic example

Here we show the input file required for the study of the solvation of a protein
by the `TMAO` solvent, which is a molecule 4 atoms. The protein is assumed to be
at infinite dilution in the simulation. The trajectory of the simulation is in `DCD`
format in this example, which is the default output of `NAMD` and `CHARMM` simulation
packages.

### Input files

The files necessary to run this would be:

- `system.pdb`: a PDB file of the complete simulated system.
- `trajectory.dcd`: the simulation trajectory, here exemplified in the DCD format.
- `script.jl`: the Julia script, described below.

These files are not provided for this example. For complete running examples, please
check our [examples](@ref examples) section.

### The Julia script

```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/basic/script.jl", String))
\`\`\`
""")
```

Given that this code is saved into a file named `script.jl`, it can be run within the Julia REPL with:
```julia
julia> include("script.jl")
```
or directly with:
```
julia -t auto script.jl
```
where `-t auto` will launch `julia` with multi-threading. It is highly recommended to use multi-threading!

!!! note
    Some newer CPUs have "fast" and "slow" cores, designed for performance or energy savings. Thus
    using all cores, with `-t auto`, may not be the best strategy for optimal performance. Experimenting
    with different number of cores using `-t N` where `N` is the number of cores used is always necessary
    for tunning performance.

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
water, a cosolvent: TMAO (trimetylamine-N-oxyde), which is a common
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
protein = select(atoms, "protein")
```
And, finally, let us use the `AtomSelection` function to setup the
structure required by the MDDF calculation:
```julia
solute = AtomSelection(protein, nmols=1)
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
tmao = select(atoms, "resname TMAO")
solvent = AtomSelection(tmao, natomspermol=14)
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
trajectory = Trajectory("trajectory.dcd", solute, solvent)
```
In the case, the trajectory is of NAMD "`DCD`" format. All formats
supported by [Chemfiles](http://chemfiles.org/Chemfiles.jl/latest/) 
are automatically recognized. 

### Finally, run the computation and get the results:

If default options are used (as the bin size of the histograms, read all
frames without skipping any), just run the `mddf` with:
```julia
results = mddf(trajectory, Options(bulk_range=(8.0, 12.0)))

```
Some optional parameters for the computation are available in the
[Options](@ref options) section. Depending on the number of atoms
and trajectory length, this can take a while. Computing a MDDF
is much more expensive than computing a regular radial distribution
function, because the normalization requires the generation of an
ideal distribution of the molecules in the system.  

### The `results` data structure obtained

The `results` data structure contains all the results of the MDDF
calculation, including:

`results.d` : Vector containing the distances to the solute. 

`results.mddf` : Vector containing the minimum-distance distribution
function at each distance.

That means, for example, that 
```julia
plot(results.d, results.mddf, xlabel="d / Å", ylabel="mddf(d)") 

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
plot(results.d, results.kb, xlabel="d / Å", ylabel="KB(d) / L / mol") 
```
to obtain:

```@raw html
<center>
<img src="../figures/kb.png" width="60%">
</center>
```

See the [Atomic and group contributions](@ref contributions) section for a
detailed account on how to obtain a molecular picture of the solvation
by splitting the MDDF in the contributions of each type of atom of the
solvent, each type of residue of the protein, etc.

### Save the results

The results can be saved into a file (with JSON format) with:
```julia
save(results, "./results.json")
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








