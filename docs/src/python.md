# [From Python](@id python)

!!! note
    Most features of the package are available through this Python interface. However, some flexibility may be reduced and, also, the tunning of the plot appearance is left to the user, as it is expected that he/she is fluent with some tools within Python if chosing this interface.

    Python 3 or greater is required.

    Please report issues, incompatibilities, or any other difficulty in using the package and its interface.
    
The following examples consider a system composed a protein solvated by a mixture of water and glycerol, built with [Packmol](http://m3g.iqm.unicamp.br/packmol). The simulations were performed with [NAMD](https://www.ks.uiuc.edu/Research/namd/) with periodic boundary conditions and a NPT ensemble at room temperature and pressure. Molecular pictures were produced with [VMD](https://www.ks.uiuc.edu/Research/vmd/).

```@raw html
<center>
<img width=50% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Data/system.png">
</center>
```

Image of the system of the example: a protein solvated by a mixture of glycreol (green) and water, at a concentration of 50%vv. The complete
example is available at [this repository](https://github.com/m3g/ComplexMixturesExamples/tree/main/Protein_in_Glycerol).

## Loading the `ComplexMixtures.py` file

The Python interface of `ComplexMixtures` is implemented in the [`ComplexMixtures.py`](./assets/ComplexMixtures.py) file. 
Just download it from the link and save it in a known path.

## Installing `juliacall`

[`juliacall`](https://github.com/cjdoris/PythonCall.jl) is a package that allows calling Julia programs from Python. Install it with

```bash
pip install juliacall
```

## Installing Julia and underlying packages

Once `juliacall` is installed, from within Python, execute:
```python
import ComplexMixtures
```
here we assume that the `ComplexMixtures.py` file is in the same directory where you launched Python.

!!! note 
     **On the first** time you execute this command, the Julia executable and the required Julia packages (`ComplexMixtures` and `PDBTools`) will be downloaded and installed. At the end of the process quit Python (not really required, but we prefer to separate the installation from the use of the module). 

## How to run this example

The [Data](https://github.com/m3g/ComplexMixturesExamples/tree/main/Protein_in_Glycerol/Data) directory contains the a pdb file of the system (`system.pdb`) and a sample from the trajectory (`glyc50.dcd`), with a few frames. It also contains the result of running the `mddf` calculation on the complete trajectory, `results_glyc50.json`. This last file was produced by `ComplexMixtures`, as indicated in the following examples. 

The sample trajectory is provided so that the first example can be run, yet do not expect that the results are the same, as the sampling is much lower in this case. The complete trajectory can be retrieved from [this link](https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing) (3GB file). 

We assume that you navigated to the directory of the example, and copied the Python module file to it: 
```bash
git clone https://github.com/m3g/ComplexMixturesExamples
cd ComplexMixturesExamples/Protein_in_Glycerol/MDDF
cp /path/to/ComplexMixtures.py ./
export JULIA_NUM_THREADS=8
```
The last line will allow Julia to execute multi-threaded, which will improve a lot the performance on most machines. Set the number of threads to the number of cores of your computer.

## Minimum-Distance Distribuion function

Note that the example here follows an identical syntax to the Julia example, except that we qualify the name of the loaded module and implicitly load the `PDBTools` package.

The script to compute the MDDFs as associated data from within python is, then:

```python
import ComplexMixtures as cm

# Load the pdb file of the system using `PDBTools`:
atoms = cm.readPDB("../Data/system.pdb")

# Create arrays of atoms with the protein and Glycerol atoms, 
# using the `select` function of the `PDBTools` package:
protein = cm.select(atoms,"protein")
glyc = cm.select(atoms,"resname GLYC")

# Setup solute and solvent structures, required for computing the MDDF, 
# with `Selection` function of the `ComplexMixtures` package:
solute = cm.Selection(protein,nmols=1)
solvent = cm.Selection(glyc,natomspermol=14)

# Read and setup the Trajectory structure required for the computations:
trajectory = cm.Trajectory("../Data/glyc50_complete.dcd",solute,solvent)

# Run the calculation and get results:
results = cm.mddf(trajectory)

# Save the reults to recover them later if required
cm.save(results,"./glyc50.json")
```

!!! note
    To change the options of the calculation, set the `Options` structure accordingly and pass it as a parameter to `mddf`. For example:
    ```python
    options = cm.Options(cutoff=10.)
    results = cm.mddf(trajectory,options)
    ```
    The complete set of options available is described [here](@ref options).


The trajectory that was loaded was for a toy-example. The complete trajectory is available [here](https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing), but it is a 3GB file. The same procedure above was performed with that file and produced the `results_Glyc50.json` file, which is available in the Data directory here. We will continue with this file instead. 

### Produce plots

#### MDDF and Kirkwood-Buff integrals

```python
import ComplexMixtures as cm
import matplotlib.pyplot as plt

# Load the actual results obtained with the complete simulation:
results = cm.load("../Data/results_glyc50.json")

# Plot MDDF and KB
fig, axs = plt.subplots(2)
axs[0].plot(results.d, results.mddf)
axs[0].set(ylabel="MDDF")

# Plot KB integral
axs[1].plot(results.d, results.kb)
axs[1].set(xlabel="distance / Angs", ylabel="MDDF")

plt.savefig("mddf_kb.png")
```

### Atomic contributions to the MDDF

Selecting the atoms corresponding to the hydroxyl groups, and of the aliphatic carbons of Glycerol. Here we list the types of the atoms as specified by the force-field. 

```python
import ComplexMixtures as cm
import matplotlib.pyplot as plt

atoms = cm.readPDB("../Data/system.pdb")
protein = cm.select(atoms,"protein")
glyc = cm.select(atoms,"resname GLYC")
solute = cm.Selection(protein,nmols=1)
solvent = cm.Selection(glyc,natomspermol=14)

# load results
results = cm.load("../Data/results_glyc50.json")

# Select atoms by name
hydroxyls = cm.list(["O1","O2","O3","H1","H2","H3"])
aliphatic = cm.list(["C1","C2","HA","HB","HC","HD"])

# Extract the contributions of the groups above
hydr_contrib = cm.contrib(solvent,results.solvent_atom,hydroxyls)
aliph_contrib = cm.contrib(solvent,results.solvent_atom,aliphatic)

# Plot
plt.plot(results.d, results.mddf)
plt.plot(results.d, hydr_contrib)
plt.plot(results.d, aliph_contrib)
plt.xlabel("distance / Angs")
plt.ylabel("MDDF")
plt.savefig("group_contributions.png")
```

!!! note
    The syntax here diverges from the Julia-only examples by requiring the lists of names
    to be converted to Julia arrays, which happens by using the `cm.list(python_list)` function calls.

                                                                                                      

































