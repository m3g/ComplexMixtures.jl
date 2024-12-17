@testitem "show methods" begin

    function test_show(
        x, s::String; 
        f64 = (x1,x2) -> isapprox(x1,x2,rtol=1e-3),
        i64 = (x1,x2) -> x1 == x2, 
        rep = Dict{String,String}(),
        vector_simplify = true,
    )
        match(f,x1,x2) = begin
            if !f(x1,x2)
                println("show method test failed with $x1 ($(typeof(x1))) == $x2 ($(typeof(x2)))")
                return false
            end
            return true
        end
        buff = IOBuffer()
        show(buff, MIME"text/plain"(), x)
        ss = String(take!(buff))
        s = replace(s, rep...)
        ss = replace(ss, rep...)
        if vector_simplify
            s = replace(s, r"\[(\S+).* (\S+)\]" => s"[ \1 \2 ]")
            ss = replace(ss, r"\[(\S+).* (\S+)\]" => s"[ \1 \2 ]")
        end
        sfields = split(s)
        ssfields = split(ss)
        all_match = true
        for (f1, f2) in zip(sfields, ssfields)
            !all_match && break
            if ispath(f1)
                all_match = last(split(f1)) == last(split(f2))
                continue
            end
            value = tryparse(Int, f1)
            if !isnothing(value)
                all_match = match(i64, value, tryparse(Int, f2))
                continue
            end
            value = tryparse(Float64, f1)
            if !isnothing(value)
                all_match = match(f64, value, tryparse(Float64,f2))
                continue
            end
            all_match = match(isequal, f1, f2)
        end
        return all_match
    end

    using ComplexMixtures
    using PDBTools: readPDB, select
    using ComplexMixtures.Testing: data_dir

    # Test simple three-molecule system: cross correlation
    atoms = readPDB("$data_dir/toy/cross.pdb")
    protein = AtomSelection(select(atoms, "protein and model 1"), nmols=1)
    water = AtomSelection(select(atoms, "resname WAT and model 1"), natomspermol=3)
    trajectory_file = "$data_dir/toy/cross.pdb"
    trajectory_format = "PDBTraj"

    @test test_show(
        Options(),
        """
        --------------------------------------------------------------------------------
        Options - ComplexMixtures 
        --------------------------------------------------------------------------------
        Trajectory frames: 
            First frame to be considered: firstframe = 1 
            Last frame to be considered (-1 is last): lastframe = -1
            Stride: stride = 1
       
        Bulk region, cutoff, and histogram:
            Bin step of histogram: binstep = 0.02
            Bulk range: >= 10.0
            (dbulk = 10.0, cutoff = 10.0, usecutoff = false)
       
        Computation details: 
            Reference atom for random rotations: irefatom = -1
            Number of random samples per frame: n_random_samples = 10
            Linked cell partition: lcell = 1
            Force garbage collection: GC = true
            Memory threshold for GC: GC_threshold = 0.3
            Seed for random number generator: 321
            Use stable random number generator: StableRNG = false 
            Number of threads to use (0 is all): nthreads = 0
            Silent output: false
        --------------------------------------------------------------------------------
        """
    )

    @test test_show(
        protein,
        """
        AtomSelection 
            1 atoms belonging to 1 molecule(s).
            Atoms per molecule: 1
            Number of groups: 1
        """
    )

    nthreads = 1
    lastframe = 2
    low_memory = true

    options = Options(;
        seed=321,
        StableRNG=true,
        nthreads,
        silent=true,
        n_random_samples=10^5,
        lastframe,
    )
    R = coordination_number(trajectory_file, protein, water, options; trajectory_format, low_memory)

    @test test_show(
        R,
        """
        --------------------------------------------------------------------------------
        MDDF Overview - ComplexMixtures - Version 2.11.4-DEV
        --------------------------------------------------------------------------------
        
        Solvent properties:
        -------------------
        
        Simulation concentration: 0.18450433782524045 mol L⁻¹
        Molar volume: 5419.9267713 cm³ mol⁻¹
        
        Concentration in bulk: 0.0 mol L⁻¹
        Molar volume in bulk: Inf cm³ mol⁻¹
        
        Solute properties:
        ------------------
        
        Simulation Concentration: 0.061501445941746814 mol L⁻¹
        Estimated solute partial molar volume: -Inf cm³ mol⁻¹
        
        Bulk range: >= 10.0 Å
        Molar volume of the solute domain: 0.0 cm³ mol⁻¹
        
        Auto-correlation: false
        
        Trajectory files and weights:
        
           /home/leandro/.julia/dev/ComplexMixtures/test/data/toy/cross.pdb - w = 1.0
        
        Long range MDDF mean (expected 1.0): 0.0 ± 0.0
        Long range RDF mean (expected 1.0): 0.0 ± 0.0
        
        --------------------------------------------------------------------------------
        """
    )

    @test test_show(
        R.volume,
        """
        Total volume: 27000.0
        Bulk volume: 0.0
        Domain volume: 0.0
        Shell volumes: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        """
    )

    @test test_show(
        R.density,
        """
        Density of solute: 3.7037037037037037e-5
        Density of solvent: 0.00011111111111111112
        Density of solvent in bulk: 0.0 
        """
    )

    @test test_show(
        SoluteGroup(select(atoms, "protein and residue 2")),
        """
        SoluteGroup defined by:
        atom_indices: [ 10 ] - 1 atoms
        """
    )

    @test test_show(
        SolventGroup(select(atoms, "protein and residue 2")),
        """
        SolventGroup defined by:
        atom_indices: [ 10 ] - 1 atoms
        """
    )

    @test test_show(
        Trajectory(trajectory_file, protein, water),
        """
        Trajectory in PDB format with:
            2 frames.
            Solute contains 1 atoms.
            Solvent contains 9 atoms.
            Unit cell in current frame: [  30.00 0 0; 0  30.00 0; 0 0  30.00 ]
        """
    )



end