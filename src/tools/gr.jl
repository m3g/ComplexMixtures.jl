"""
    gr(r::AbstractVector{<:Real}, count::AbstractVector{<:Real}, density::Real, binstep::Real)

Computes the radial distribution function from the count data and the density.

This is exactly a conventional g(r) if a single atom was chosen as the solute and solvent selections.

Returns both the g(r) and the kb(r)

"""
function gr(r::AbstractVector{<:Real}, count::AbstractVector{<:Real}, density::Real, binstep::Real)
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
    gr(R::Result)

If a Result structure is provided without further details, use the rdf count and the bulk solvent density.

If the solute contains a single atom per molecule, the distances in `R.rdf_count` are radial
distances to the reference atom of the solvent molecules, and the g(r) and KB integral are
computed with the analytical spherical shell volumes.

If the solute contains more than one atom per molecule, the distances are minimum distances
from the solute to the reference atom of the solvent, and the volume elements are not spherical
shells. In this case, the distribution is normalized by the random reference state (`R.rdf`), 
and the KB integral is `R.kb_rdf`.

"""
function gr(R::Result)
    if R.solute.natomspermol == 1
        return gr(R.d, R.rdf_count, R.density.solvent_bulk, R.files[1].options.binstep)
    else
        return copy(R.rdf), copy(R.kb_rdf)
    end
end

@testitem "Radial distribution" begin
    using ComplexMixtures: gr, mddf, Trajectory, Options, AtomSelection
    using ComplexMixtures: data_dir
    using PDBTools: read_pdb, select
    atoms = read_pdb("$data_dir/NAMD/structure.pdb")
    options = Options(seed=321, StableRNG=true, nthreads=1, silent=true)
    OH2 = AtomSelection(select(atoms, "water and name OH2"), natomspermol=1)
    traj = Trajectory("$data_dir/Gromacs/trajectory.xtc", OH2)
    R = mddf(traj, options)
    gr1, kb1 = gr(R)
    @test R.rdf_count ≈ R.md_count
    @test gr1[end] ≈ 1.0 rtol = 0.1
    @test kb1[end] ≈ 20.0 rtol = 0.1
end

@testitem "Radial distribution - multi-atom solute" begin
    using ComplexMixtures: gr, load
    using ComplexMixtures: data_dir
    R = load("$data_dir/NAMD/protein_tmao.json")
    @test R.solute.natomspermol > 1
    g, kb = gr(R)
    @test g == R.rdf
    @test kb == R.kb_rdf
    @test g !== R.rdf
end
