# [Save and load results](@id save)

Three functions serve the purpose of saving and loading the results
obtained with MDDF:

## Save data to recover it later 

```julia
MDDF.save(results,"results.json")
```
where `results` is the output data structure of the `MDDF.mddf()`
calculation, and `results.json` is the output file to be created. The
file is written in `JSON` format, thus is not naturally human-readable.

## Recover saved data

```julia
results = MDDF.read("results.json")
```
The `MDDF.read` function reads the output of the `save` function above,
and restores the results data structure.

## Write data in a human-readable format

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







