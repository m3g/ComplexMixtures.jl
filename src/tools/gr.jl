"""
    gr(r::Vector{Float64}, count::Vector{Float64}, density::Float64, binstep::Float64)

Computes the radial distribution function from the count data and
the density.

This is exactly a conventional g(r) 
if a single atom was chosen as the solute and solvent selections.

Returns both the g(r) and the kb(r)

"""
function gr(r::Vector{Float64}, count::Vector{Float64}, density::Float64, binstep::Float64)
    nbins = length(r)
    gr = zeros(nbins)
    kb = zeros(nbins)
    for i = 1:nbins
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
    gr(R::Result) = gr(R.d,R.rdf_count,R.density.solvent_bulk,R.options.binstep)

If a Result structure is provided without further details, use the rdf count
and the bulk solvent density.
"""
gr(R::Result) = gr(R.d, R.rdf_count, R.density.solvent_bulk, R.options.binstep)

@testitem "Radial distribution" begin
    using ComplexMixtures, PDBTools
    using ComplexMixtures.Testing
    const CM = ComplexMixtures
    dir = "$(Testing.data_dir)/Gromacs"
    atoms = readPDB("$dir/system.pdb")
    options = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true)
    emi = AtomSelection(select(atoms, "resname EMI and name H01"), natomspermol = 20)
    traj = Trajectory("$dir/trajectory.xtc", emi)
    R = mddf(traj, options)
    gr1, kb1 = gr(R)
    @test isapprox(gr1, R.mddf)
    @test isapprox(kb1, R.kb)
end
