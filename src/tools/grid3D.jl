"""
    grid3D(
        result::Result, atoms, output_file::Union{Nothing,String} = nothing; 
        dmin=1.5, dmax=5.0, step=0.5, silent = false, type = :mddf,
    )

This function builds the grid of the 3D density function and fills an array of
mutable structures of type Atom, containing the position of the atoms of 
grid, the closest atom to that position, and distance. 

## Positional arguments

- `result` is a `ComplexMixtures.Result` object 
- `atoms` is a vector of `PDBTools.Atom`s with all the atoms of the system. 
- `output_file` is the name of the file where the grid will be written. If `nothing`, the grid is only returned as a matrix. 

## Keyword (optional) arguments

- `dmin` and `dmax` define the range of distance where the density grid will be built, and `step`
    defines how fine the grid must be. Be aware that fine grids involve usually a very large (hundreds
    of thousands points).
- `silent` is a boolean to suppress the progress bar.
- `type` can be `:mddf`, `:coordination_number`, or `:md_count`, depending on the data available or desired from the results.

The output PDB has the following characteristics:

- The positions of the atoms are grid points. 
- The identity of the atoms correspond to the identity of the protein atom contributing to the property at that point (the closest protein atom). 
- The temperature-factor column (`beta`) contains the relative contribution of that atom to the property at the corresponding distance. 
- The `occupancy` field contains the distance itself.

### Example

```julia-repl
julia> using ComplexMixtures, PDBTools

julia> atoms = read_pdb("./system.pdb");

julia> R = ComplexMixtures.load("./results.json");

julia> grid = grid3D(R, atoms, "grid.pdb");
```

Examples of how the grid can be visualized are provided in the user guide of `ComplexMixtures`. 

"""
function grid3D(
    result::Result,
    atoms,
    output_file::Union{Nothing,String}=nothing;
    dmin=1.5,
    dmax=5.0,
    step=0.5,
    silent=false,
    type=:mddf,
)

    if result.solute.custom_groups
        throw(ArgumentError("""\n

            The 3D grid can only be built if the contributions of all atoms of the solute are recorded.
            It is not compatible with predefined groups of atoms.
                
        """))
    end

    # Simple function to interpolate data
    interpolate(x₁, x₂, y₁, y₂, xₙ) = y₁ + (y₂ - y₁) / (x₂ - x₁) * (xₙ - x₁)

    # Maximum and minimum coordinates of the solute
    solute_atoms = atoms[result.solute.indices]
    lims = PDBTools.maxmin(solute_atoms)
    n = @. ceil(Int, (lims.xlength + 2 * dmax) / step + 1)

    # Building the grid with the nearest solute atom information
    igrid = 0
    AtomType = typeof(PDBTools.Atom()) # to support PDBTools < 2 (which does not have a Atom{T} constructor)
    grid = AtomType[]
    grid_lock = ReentrantLock()
    p = Progress(prod(n); desc="Building grid...", enabled=!silent)
    Threads.@threads for ix_inds in ChunkSplitters.chunks(1:n[1]; n=Threads.nthreads())
        for ix in ix_inds, iy in 1:n[2], iz in 1:n[3]
            next!(p)
            x = lims.xmin[1] - dmax + step * (ix - 1)
            y = lims.xmin[2] - dmax + step * (iy - 1)
            z = lims.xmin[3] - dmax + step * (iz - 1)
            rgrid = -1
            _, iat, r = PDBTools.closest(SVector(x, y, z), solute_atoms)
            if (dmin < r < dmax)
                if rgrid < 0 || r < rgrid
                    at = solute_atoms[iat]
                    # Get contribution of this atom to the MDDF
                    c = contributions(result, SoluteGroup(SVector(PDBTools.index(at),)); type)
                    # Interpolate c at the current distance
                    iright = findfirst(d -> d > r, result.d)
                    ileft = iright - 1
                    cᵣ = interpolate(
                        result.d[ileft],
                        result.d[iright],
                        c[ileft],
                        c[iright],
                        r,
                    )
                    if cᵣ > 0
                        gridpoint = AtomType(
                            index=PDBTools.index(at),
                            index_pdb=PDBTools.index_pdb(at),
                            name=PDBTools.name(at),
                            chain=PDBTools.chain(at),
                            resname=PDBTools.resname(at),
                            resnum=PDBTools.resnum(at),
                            x=x,
                            y=y,
                            z=z,
                            occup=r,
                            beta=cᵣ,
                            model=PDBTools.model(at),
                            segname=PDBTools.segname(at),
                        )
                        if rgrid < 0
                            @lock grid_lock begin
                                igrid += 1
                                push!(grid, gridpoint)
                            end
                        elseif r < rgrid
                            @lock grid_lock begin
                                grid[igrid] = gridpoint
                            end
                        end
                        rgrid = r
                    end # cᵣ>0
                end # rgrid
            end # dmin/dmax
        end # ix, iy, iz
    end # chunks

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
        PDBTools.writePDB(grid, output_file)
        silent || println("Grid written to $output_file")
    end
    return grid
end

@testitem "grid3D" begin
    using PDBTools
    using ComplexMixtures
    using ComplexMixtures: data_dir
    dir = "$data_dir/NAMD"
    atoms = read_pdb("$dir/structure.pdb")

    # Test argument error: no custom groups can be defined
    protein = AtomSelection(select(atoms, "protein"); group_atom_indices=[findall(sel"resname ARG", atoms)], nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    options = Options(
        stride=5,
        seed=321,
        StableRNG=true,
        nthreads=1,
        silent=true,
        n_random_samples=100,
    )
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
    R = mddf(traj, options)
    @test_throws ArgumentError grid3D(R, atoms, tempname())

    # Test properties of the grid around a specific residue
    solute = AtomSelection(select(atoms, "protein and residue 46"), nmols=1)
    solvent = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$dir/trajectory.dcd", solute, solvent)
    grid_file = tempname() * ".pdb"
    options = Options(
        stride=5,
        seed=321,
        StableRNG=true,
        nthreads=1,
        silent=true,
    )
    R = mddf(traj, options)
    grid = grid3D(R, atoms, grid_file)
    @test length(grid) ≈ 1539 atol = 3
    c05 = filter(at -> beta(at) > 0.5, grid)
    @test length(c05) == 14
    @test all(at -> element(at) == "O", c05)
    @test all(at -> occup(at) < 2.0, c05)

    # Test if the file was properly written
    grid_read = read_pdb(grid_file)
    for property in [:name, :resname, :chain, :resnum]
        @test all(p -> getproperty(first(p), property) == getproperty(last(p), property), zip(grid, grid_read))
    end
    for property in [:x, :y, :z, :occup, :beta]
        @test all(p -> isapprox(getproperty(first(p), property), getproperty(last(p), property), atol=1e-2), zip(grid, grid_read))
    end
    rm(grid_file)

    # Test grid generation with coordination number only
    R = coordination_number(traj, options)
    grid = grid3D(R, atoms, grid_file; type=:coordination_number)
    grid_read = read_pdb(grid_file)
    for property in [:name, :resname, :chain, :resnum]
        @test all(p -> getproperty(first(p), property) == getproperty(last(p), property), zip(grid, grid_read))
    end
    for property in [:x, :y, :z, :occup, :beta]
        @test all(p -> isapprox(getproperty(first(p), property), getproperty(last(p), property), atol=1e-2), zip(grid, grid_read))
    end
    rm(grid_file)
end

