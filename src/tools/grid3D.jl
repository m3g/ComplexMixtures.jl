"""

```
grid3D(solute,solute_atoms,mddf_result,output_file; dmin=1.5, ddax=5.0, step=0.5)
```

This function builds the grid of the 3D density function and fills an array of
mutable structures of type Atom, containing the position of the atoms of 
grid, the closest atom to that position, and distance. 

`solute` is a `ComplexMixtuers.Selection`, defining the solute. `solute_atoms` is the corresponding
vector of `PDBTools.Atom`s, and `mddf_result` is the result of a `mddf_result` calculation with 
the correspondign solute. 

`dmin` and `dmax` define the range of distance where the density grid will be built, and `step`
defines how fine the grid must be. Be aware that fine grids involve usually a very large (hundreds
of thousands points).

All parameters can be provides as keyword parameters.

### Example

```julia-repl
julia> using ComplexMixtures, PDBTools

julia> pdb = readPDB("./system.pdb");

julia> R = ComplexMixtures.load("./results.json");

julia> protein = select(pdb,"protein");

julia> solute = ComplexMixtures.Selection(protein,nmols=1);

julia> grid = ComplexMixtures.grid3D(solute=solute, solute_atoms=protein, mddf_result=R, output_file="grid.pdb");

```

`grid` will contain a vector of `Atom`s with the information of the MDDF at each grid point, and the
same data will be written in the `grid.pdb` file. This PDB file can be opened in VMD, for example, and contain
in the `beta` field the contribution of each protein residue to the MDDF at each point in space relative to the 
protein, and in the `occupancy` field the distance to the protein. Examples of how this information can be
visualized are provided in the user guide of `ComplexMixtures`. 


"""
grid3D(;
    solute = nothing,
    solute_atoms = nothing,
    mddf_result = nothing,
    output_file = nothing,
    dmin = 1.5,
    dmax = 5.0,
    step = 0.5,
) = grid3D(
    solute,
    solute_atoms,
    mddf_result,
    output_file,
    dmin = dmin,
    dmax = dmax,
    step = step,
)

function grid3D(
    solute,
    solute_atoms,
    mddf_result,
    output_file;
    dmin = 1.5,
    dmax = 5.0,
    step = 0.5,
)

    if nothing in (solute, solute_atoms, mddf_result)
        error("grid3D requires `solute`, `solute_atoms` and `mddf_result` definitions.")
    end

    # Simple function to interpolate data
    interpolate(x₁, x₂, y₁, y₂, xₙ) = y₁ + (y₂ - y₁) / (x₂ - x₁) * (xₙ - x₁)

    # Maximum and minimum coordinates of the solute
    lims = maxmin(solute_atoms)
    n = @. ceil(Int, (lims.xlength + 2 * dmax) / step + 1)

    # Building the grid with the nearest solute atom information
    igrid = 0
    grid = PDBTools.Atom[]
    for ix = 1:n[1], iy = 1:n[2], iz = 1:n[3]
        x = lims.xmin[1] - dmax + step * (ix - 1)
        y = lims.xmin[2] - dmax + step * (iy - 1)
        z = lims.xmin[3] - dmax + step * (iz - 1)
        rgrid = -1
        _, iat, r = PDBTools.closest(SVector(x, y, z), solute_atoms)
        if (dmin < r < dmax)
            if rgrid < 0 || r < rgrid
                at = solute_atoms[iat]
                # Get contribution of this atom to the MDDF
                c = ComplexMixtures.contrib(
                    solute,
                    mddf_result.solute_atom,
                    [at.index_pdb],
                )
                # Interpolate c at the current distance
                iright = findfirst(d -> d > r, mddf_result.d)
                ileft = iright - 1
                cᵣ = interpolate(
                    mddf_result.d[ileft],
                    mddf_result.d[iright],
                    c[ileft],
                    c[iright],
                    r,
                )
                if cᵣ > 0
                    gridpoint = Atom(
                        index = at.index,
                        index_pdb = at.index_pdb,
                        name = at.name,
                        chain = at.chain,
                        resname = at.resname,
                        resnum = at.resnum,
                        x = x,
                        y = y,
                        z = z,
                        occup = r,
                        beta = cᵣ,
                        model = at.model,
                        segname = at.segname,
                    )
                    if rgrid < 0
                        igrid += 1
                        push!(grid, gridpoint)
                    elseif r < rgrid
                        grid[igrid] = gridpoint
                    end
                    rgrid = r
                end # cᵣ>0
            end # rgrid
        end # dmin/dmax
    end #ix

    # Now will scale the density to be between 0 and 99.9 in the temperature
    # factor column, such that visualization is good enough
    bmin, bmax = +Inf, -Inf
    for gridpoint in grid
        bmin = min(bmin, gridpoint.beta)
        bmax = max(bmax, gridpoint.beta)
    end
    for gridpoint in grid
        gridpoint.beta = (gridpoint.beta - bmin) / (bmax - bmin)
    end

    if !isnothing(output_file)
        writePDB(grid, output_file)
    end
    return grid
end
