"""
    grid3D(
        result::Result, atoms, output_file::Union{Nothing,String} = nothing; 
        dmin=1.5, ddax=5.0, step=0.5, silent = false,
    )

This function builds the grid of the 3D density function and fills an array of
mutable structures of type Atom, containing the position of the atoms of 
grid, the closest atom to that position, and distance. 

`result` is a `ComplexMixtures.Result` object 
`atoms` is a vector of `PDBTools.Atom`s with all the atoms of the system. 
`output_file` is the name of the file where the grid will be written. If `nothing`, the grid is only returned as a matrix. 

`dmin` and `dmax` define the range of distance where the density grid will be built, and `step`
defines how fine the grid must be. Be aware that fine grids involve usually a very large (hundreds
of thousands points).

`silent` is a boolean to suppress the progress bar.

### Example

```julia-repl
julia> using ComplexMixtures, PDBTools

julia> atoms = readPDB("./system.pdb");

julia> R = ComplexMixtures.load("./results.json");

julia> grid = grid3D(R, atoms, "grid.pdb");
```

`grid` will contain a vector of `Atom`s with the information of the MDDF at each grid point, and the
same data will be written in the `grid.pdb` file. This PDB file can be opened in VMD, for example, and contain
in the `beta` field the contribution of each protein residue to the MDDF at each point in space relative to the 
protein, and in the `occupancy` field the distance to the protein. Examples of how this information can be
visualized are provided in the user guide of `ComplexMixtures`. 

"""
function grid3D(
    result::Result, 
    atoms, 
    output_file::Union{Nothing,String} = nothing; 
    dmin=1.5, 
    dmax=5.0, 
    step=0.5,
    silent = false,
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
    grid = PDBTools.Atom[]
    grid_lock = ReentrantLock()
    silent || (p = Progress(prod(n), "Building grid..."))
    Threads.@threads for ix_inds in ChunkSplitters.chunks(1:n[1]; n=Threads.nthreads())
        for ix in ix_inds, iy in 1:n[2], iz in 1:n[3]
            silent || next!(p)
            x = lims.xmin[1] - dmax + step * (ix - 1)
            y = lims.xmin[2] - dmax + step * (iy - 1)
            z = lims.xmin[3] - dmax + step * (iz - 1)
            rgrid = -1
            _, iat, r = PDBTools.closest(SVector(x, y, z), solute_atoms)
            if (dmin < r < dmax)
                if rgrid < 0 || r < rgrid
                    at = solute_atoms[iat]
                    # Get contribution of this atom to the MDDF
                    c = contributions(result, SoluteGroup(SVector(PDBTools.index(at),)))
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
                        gridpoint = PDBTools.Atom(
                            index = PDBTools.index(at),
                            index_pdb = PDBTools.index_pdb(at),
                            name = PDBTools.name(at),
                            chain = PDBTools.chain(at),
                            resname = PDBTools.resname(at),
                            resnum = PDBTools.resnum(at),
                            x = x,
                            y = y,
                            z = z,
                            occup = r,
                            beta = cᵣ,
                            model = PDBTools.model(at),
                            segname = PDBTools.segname(at),
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
    using ComplexMixtures.Testing: data_dir
    dir = "$data_dir/NAMD"
    atoms = readPDB("$dir/structure.pdb")

    # Test argument error: no custom groups can be defined
    protein = AtomSelection(select(atoms, "protein"); group_atom_indices = [ findall(sel"resname ARG", atoms) ], nmols = 1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)
    options = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true)
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
    R = mddf(traj, options)
    @test_throws ArgumentError grid3D(R, atoms, tempname())

    # Test properties of the grid around a specific residue
    solute = AtomSelection(select(atoms, "protein and residue 46"), nmols = 1)
    solvent = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$dir/trajectory.dcd", solute, solvent)
    grid_file = tempname()*".pdb"
    R = mddf(traj, options)
    grid = grid3D(R, atoms, grid_file)
    @test length(grid) == 1539
    c05 = filter(at -> beta(at) > 0.5, grid)
    @test length(c05) == 15
    @test all(at -> element(at) == "O", c05)
    @test all(at -> occup(at) < 1.8, c05)

    # Test if the file was properly written
    grid_read = readPDB(grid_file)
    for property in [:name, :resname, :chain, :resnum ]
        @test all(p -> getproperty(first(p), property) == getproperty(last(p), property), zip(grid, grid_read))
    end
    for property in [ :x, :y, :z, :occup, :beta ]
        @test all(p -> isapprox(getproperty(first(p), property),getproperty(last(p), property),atol=1e-2), zip(grid, grid_read))
    end
end

