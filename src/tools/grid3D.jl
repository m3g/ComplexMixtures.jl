"""
    grid3D(
        result::Result, atoms, output_file::Union{Nothing,String} = nothing; 
        dmin=1.5, dmax=5.0, step=0.5, silent = false, type = :mddf, molecule = 1,
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
- `molecule`: if the solute is composed of more than one molecule (for example, glycerol in a water/glycerol mixture), 
    the grid is built around this molecule of the solute (the first one, by default). The contributions are those of 
    each atom of a single solute molecule, averaged over all molecules, and the grid is built around the structure of
    the chosen molecule in `atoms`.

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
    molecule::Integer=1,
)

    if result.solute.custom_groups
        throw(ArgumentError("""\n

            The 3D grid can only be built if the contributions of all atoms of the solute are recorded.
            It is not compatible with predefined groups of atoms.
                
        """))
    end

    # Simple function to interpolate data
    interpolate(x₁, x₂, y₁, y₂, xₙ) = y₁ + (y₂ - y₁) / (x₂ - x₁) * (xₙ - x₁)

    # Atoms of the solute molecule around which the grid is built
    (; nmols, natomspermol, indices) = result.solute
    1 <= molecule <= nmols || throw(ArgumentError("molecule must be between 1 and the number of solute molecules ($nmols). Got: $molecule"))
    solute_atoms = atoms[indices[(molecule-1)*natomspermol+1:molecule*natomspermol]]

    # Maximum and minimum coordinates of the solute
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

"""
    volumetric_data(grid::AbstractVector{<:PDBTools.Atom}; sigma=0.0, step=nothing, value=PDBTools.beta)

Converts the grid computed by `grid3D` into volumetric data (a `PDBTools.VolumetricData` object), that is,
the values of the contributions on a regular (dense) grid. The positions of the regular grid that are not
grid points of `grid3D` are filled with zeros. 

The volumetric data can be displayed as isosurfaces with `PDBTools.visualize`, or written to a file in the 
OpenDX (`.dx`) format with `PDBTools.write_dx`, to be loaded in other visualization software 
(VMD, PyMOL, ChimeraX).

## Keyword (optional) arguments

- `sigma`: width (in Å) of the Gaussian function used to smooth the data. With `sigma=0` (default), the data is not smoothed.
  The values of the grid of `grid3D` are the contributions of the closest solute atom to each point, thus they change
  abruptly between neighboring points that are closest to different atoms. Smoothing (with `sigma` of the order 
  of the grid `step`) produces continuous isosurfaces. Smoothing preserves the sum of the values, thus the maximum
  value decreases with `sigma`.
- `step`: the spacing of the grid. By default, it is obtained from the positions of the grid points, which must be 
  the same `step` used in `grid3D`.
- `value`: the function that returns the value associated to each grid point. By default, it is the
  `beta` field, which contains the relative contribution of the closest solute atom to the property at that point.

### Example

```julia-repl
julia> using ComplexMixtures, PDBTools

julia> atoms = read_pdb("./system.pdb");

julia> R = ComplexMixtures.load("./results.json");

julia> grid = grid3D(R, atoms; dmin=1.5, dmax=3.5);

julia> density = volumetric_data(grid; sigma=0.5);

julia> write_dx("density.dx", density) # write the data to a file

julia> visualize(select(atoms, "protein") => (;), density => (isovalue=0.05,)) # isosurface of the data
```

!!! compat
    This function is available in ComplexMixtures 2.19.0 or greater, and requires PDBTools 3.41.0 or greater.

"""
function volumetric_data(
    grid::AbstractVector{<:PDBTools.Atom};
    sigma::Real=0.0,
    step::Union{Nothing,Real}=nothing,
    value::Function=PDBTools.beta,
)
    isempty(grid) && throw(ArgumentError("The grid is empty."))
    sigma >= 0 || throw(ArgumentError("sigma must be non-negative. Got: $sigma"))
    positions = [SVector(at.x, at.y, at.z) for at in grid]
    if isnothing(step)
        # Smallest separation between the coordinates of the grid points
        step = Inf
        for k in 1:3
            x = sort!(unique(round.(getindex.(positions, k); digits=4)))
            for i in 2:length(x)
                step = min(step, x[i] - x[i-1])
            end
        end
        isfinite(step) || throw(ArgumentError("Could not determine the grid step. Provide it with the `step` keyword."))
        # The coordinates are stored in single precision
        step = round(step; digits=3)
    end
    step > 0 || throw(ArgumentError("step must be positive. Got: $step"))
    # Kernel of the Gaussian smoothing, and padding of the grid, such that the data goes to zero at the borders
    nkernel = sigma > 0 ? ceil(Int, 3 * sigma / step) : 0
    npad = nkernel + 1
    xmin = reduce((x, y) -> min.(x, y), positions) .- npad * step
    xmax = reduce((x, y) -> max.(x, y), positions) .+ npad * step
    n = round.(Int, (xmax .- xmin) ./ step) .+ 1
    data = zeros(n...)
    for (at, x) in zip(grid, positions)
        f = (x .- xmin) ./ step
        i = round.(Int, f)
        if any(abs.(f .- i) .> 1e-2)
            throw(ArgumentError("Grid point at $x is not on a regular grid with step $step."))
        end
        data[(i .+ 1)...] = value(at)
    end
    if sigma > 0
        weights = [exp(-((k * step)^2) / (2 * sigma^2)) for k in -nkernel:nkernel]
        weights ./= sum(weights)
        # Separable convolution: the Gaussian is applied in each dimension
        smoothed = similar(data)
        for dim in 1:3
            _smooth_dimension!(smoothed, data, weights, nkernel, dim)
            data, smoothed = smoothed, data
        end
    end
    return PDBTools.VolumetricData(data; origin=xmin, step=step)
end

# Convolution of the map with the weights along dimension dim
function _smooth_dimension!(smoothed::Array{T,3}, map::Array{T,3}, weights, nkernel, dim) where {T}
    unit = CartesianIndex(ntuple(d -> d == dim ? 1 : 0, 3))
    for I in CartesianIndices(map)
        s = zero(T)
        for (k, w) in zip(-nkernel:nkernel, weights)
            J = I + k * unit
            checkbounds(Bool, map, J) || continue
            s += w * map[J]
        end
        smoothed[I] = s
    end
    return smoothed
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


@testitem "volumetric_data" begin
    using PDBTools
    using ComplexMixtures
    using ComplexMixtures: data_dir
    # Points on a regular grid
    grid = [Atom(x=0.5 * i, y=0.5 * j, z=1.0, beta=float(i + j)) for i in 1:3 for j in 1:2]
    v = volumetric_data(grid)
    @test v isa PDBTools.VolumetricData
    @test size(v.data) == (5, 4, 3)
    @test v.origin ≈ [0.0, 0.0, 0.5]
    @test v.step ≈ [0.5, 0.5, 0.5]
    @test sum(v.data) ≈ sum(beta.(grid))
    # point (i=1, j=1) is at index (1,1,1) + padding of 1
    @test v.data[2, 2, 2] == 2.0
    # Smoothing preserves the sum of the values (the grid is padded)
    vs = volumetric_data(grid; sigma=0.5)
    @test sum(vs.data) ≈ sum(beta.(grid))
    @test maximum(vs.data) < maximum(beta.(grid))
    # Other properties
    @test sum(volumetric_data(grid; value=at -> 2 * beta(at)).data) ≈ 2 * sum(beta.(grid))
    # Write and read DX file
    dx_file = tempname() * ".dx"
    write_dx(dx_file, vs)
    @test read_dx(dx_file).data ≈ vs.data rtol = 1e-5
    rm(dx_file)
    # Errors
    @test_throws ArgumentError volumetric_data(PDBTools.Atom[])
    @test_throws ArgumentError volumetric_data(grid; sigma=-1.0)
    @test_throws ArgumentError volumetric_data([grid; Atom(x=0.7, y=0.5, z=1.0)]; step=0.5)
    # Grid from a simulation
    atoms = read_pdb("$data_dir/NAMD/structure.pdb")
    solute = AtomSelection(select(atoms, "protein and residue 46"), nmols=1)
    solvent = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$data_dir/NAMD/trajectory.dcd", solute, solvent)
    R = mddf(traj, Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true))
    grid = grid3D(R, atoms; silent=true)
    v = volumetric_data(grid; sigma=0.5)
    @test v.step ≈ [0.5, 0.5, 0.5]
    @test sum(v.data) ≈ sum(beta.(grid)) rtol = 1e-4
end

@testitem "grid3D with multiple solute molecules" begin
    using PDBTools
    using ComplexMixtures
    using ComplexMixtures: data_dir
    atoms = read_pdb("$data_dir/NAMD/structure.pdb")
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$data_dir/NAMD/trajectory.dcd", tmao, water)
    R = mddf(traj, Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true))
    grid1 = grid3D(R, atoms; silent=true)
    # The grid is built around a single molecule: the closest atoms belong to the first molecule
    first_molecule = atoms[R.solute.indices[1:14]]
    @test all(p -> index(p) in index.(first_molecule), grid1)
    grid2 = grid3D(R, atoms; silent=true, molecule=2)
    second_molecule = atoms[R.solute.indices[15:28]]
    @test all(p -> index(p) in index.(second_molecule), grid2)
    @test_throws ArgumentError grid3D(R, atoms; silent=true, molecule=0)
    @test_throws ArgumentError grid3D(R, atoms; silent=true, molecule=R.solute.nmols + 1)
end
