"""
    contourf_per_residue(results::Result, atoms::AbstractVector{PDBTools.Atom}; 
        residue_range=1:length(eachresidue(atoms)), 
        dmin=1.5, dmax=3.5, 
        oneletter=false,
        xlabel="Residue",
        ylabel="r / Å"
    )

Plot the contribution of each residue to the solute-solvent pair distribution function as a contour plot.
This function requires loading the `Plots` and `PDBTools` packages.

# Arguments

- `results::Result`: The result of a `mddf` call.
- `atoms::AbstractVector{Atom}`: The atoms of the solute.

# Optional arguments

- `residue_range::UnitRange{Int}`: The range of residues to plot. Default is `1:length(eachresidue(atoms))`.
- `dmin::Real`: The minimum distance to plot. Default is `1.5`.
- `dmax::Real`: The maximum distance to plot. Default is `3.5`.
- `oneletter::Bool`: Use one-letter residue codes. Default is `false`. One-letter codes are only available for the 20 standard amino acids.
- `xlabel` and `ylabel`: Labels for the x and y axes. Default is `"Residue"` and `"r / Å"`.

# Example

```julia-repl
julia> using ComplexMixtures, Plots, PDBTools

julia> results = load("mddf.json")

julia> atoms = readPDB("system.pdb", "protein")

julia> plt = contourf_per_residue(results, atoms; oneletter=true)
```

This will produce a plot with the contribution of each residue to the solute-solvent pair distribution function,
as a contour plot, with the residues in the x-axis and the distance in the y-axis.

The resulting plot can customized using the standard mutating `plot!` function, for example, 

```julia-repl
julia> plot!(plt, size=(800, 400), title="Contribution per residue")
```

!!! compat
    This function requires loading the `Plots` package and is available in 
    ComplexMixtures v2.2.0 or greater.

"""
function contourf_per_residue(args...; kwargs...)
    throw(ArgumentError("""\n

        This function requires loading the `Plots` and `PDBTools` packages.
    
    """))
end