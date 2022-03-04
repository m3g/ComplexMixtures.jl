"""

```
cutoffdistances!(cutoff::Float64,
                 x_solute::AbstractVector{T},
                 x_solvent::Vector{T},
                 lc_solvent::LinkedCells,
                 box::Box, 
                 dc::CutoffDistances) where T
```

This routine that returns a list of the distances 
between atoms that are smaller than a specified cutoff,
for a given set of coordinates.

Returns `nd`, the number of distances smaller than the cutoff, and modifies `dc`.

"""
function cutoffdistances!(
    cutoff::Float64,
    x_solute::AbstractVector{T},
    x_solvent::Vector{T},
    lc_solvent::LinkedCells,
    box::Box,
    dc::CutoffDistances,
) where {T}

    # Reset the dc structure 
    reset!(dc)

    # Loop over solute atoms
    for iat in eachindex(x_solute)
        xat = x_solute[iat]
        # Check the cell of this atom
        i, j, k = icell3D(xat, box)
        # Loop over vicinal cells to compute distances to solvent atoms, and
        # add data to dc structure (includes current cell)
        for ic = i-box.lcell:i+box.lcell
            for jc = j-box.lcell:j+box.lcell
                for kc = k-box.lcell:k+box.lcell
                    cutoffdcell!(
                        cutoff,
                        iat,
                        xat,
                        x_solvent,
                        lc_solvent,
                        box,
                        ic,
                        jc,
                        kc,
                        dc,
                    )
                end
            end
        end
    end

    return nothing
end
