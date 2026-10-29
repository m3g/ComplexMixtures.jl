#=
    _gr(r::AbstractVector{<:Real}, count::AbstractVector{<:Real}, density::Real, binstep::Real)

Computes the radial distribution function, g(r), and the corresponding Kirkwood-Buff integral
(in cm³/mol), from the number of site counts at each distance, `count`, normalized by the volumes
of the spherical shells (of width `binstep`) and by the `density` of the solvent.

Returns the tuple `(g, kb)`.

=#
function _gr(r::AbstractVector{<:Real}, count::AbstractVector{<:Real}, density::Real, binstep::Real)
    nbins = length(r)
    gr = zeros(nbins)
    kb = zeros(nbins)
    for i in eachindex(r)
        gr[i] = (count[i] / sphericalshellvolume(i, binstep)) / density
        if i == 1
            kb[i] = 4π * (gr[i] - 1) * r[i]^2 * binstep
        else
            kb[i] = kb[i-1] + 4π * (gr[i] - 1) * r[i]^2 * binstep
        end
    end
    @. kb = units.Angs3tocm3permol * kb
    return gr, kb
end


"""
    gr(R::MDDFResult)

Computes the radial distribution function, g(r), and the corresponding Kirkwood-Buff integral
(in cm³/mol), of the reference atom of the solvent molecules (`irefatom` in [`Options`](@ref)) 
relative to the solute, from the results of a [`mddf`](@ref) calculation. 

The distance counts (`R.proximal_count`) are normalized by the volumes of the spherical shells and 
by the bulk density of the solvent (`R.density.solvent_bulk`). This is only possible if the solute 
has a single atom per molecule (for example, the oxygen atom of water molecules). 

Returns the tuple `(g, kb)`.

For solutes with more than one atom per molecule an `ArgumentError` is thrown, since the volume
elements are not spherical shells. In this case, the distribution of the reference atom of the
solvent normalized by the random reference state, and its KB integral, are available in the
`R.proximal` and `R.kb_proximal` fields.

"""
function gr(R::MDDFResult)
    if R.solute.natomspermol != 1
        throw(ArgumentError("""\n
            The radial distribution function can only be computed if the solute has a single atom 
            per molecule (the solute has $(R.solute.natomspermol) atoms per molecule). 

            The proximal distribution of the reference atom of the solvent, and its KB integral,
            are available in the `R.proximal` and `R.kb_proximal` fields.

        """))
    end
    return _gr(R.d, R.proximal_count, R.density.solvent_bulk, R.files[1].options.binstep)
end

@testitem "Radial distribution" begin
    using ComplexMixtures: gr, mddf, Trajectory, Options, AtomSelection
    using LazyArtifacts
    using ComplexMixtures: test_file
    using PDBTools: read_pdb, select
    atoms = read_pdb(test_file("protein_tmao.pdb"))
    options = Options(seed=321, stable_rng=true, nthreads=1, silent=true)
    OH2 = AtomSelection(select(atoms, "water and name OH2"), natomspermol=1)
    traj = Trajectory(test_file("gromacs_trajectory.xtc"), OH2)
    R = mddf(traj, options)
    gr1, kb1 = gr(R)
    @test R.proximal_count ≈ R.md_count
    @test gr1[end] ≈ 1.0 rtol = 0.1
    @test kb1[end] ≈ 20.0 rtol = 0.1
end

@testitem "Radial distribution - multi-atom solute" begin
    using ComplexMixtures: gr, load
    using LazyArtifacts
    using ComplexMixtures: test_file
    R = load(MDDFResult, test_file("namd_protein_tmao.json"))
    @test R.solute.natomspermol > 1
    @test_throws ArgumentError gr(R)
end
