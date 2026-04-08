"""
    kbi(R::Result; correction=:first_order)

Compute the Kirkwood-Buff Integral (KBI) from the spatial correlation data 
contained in the `Result` object `R`.

This function integrates the solvent density fluctuations around the solute. 
Because simulations are performed in finite, closed periodic boxes, the raw 
integral suffers from thermodynamic finite-size depletion effects. The 
`correction` keyword allows for systematic dampening of this boundary noise.

# Arguments
- `R::Result`: A result object containing the spatial distribution data, including 
  shell volumes, solvent bulk density, and molecular counts.

# Keyword Arguments
- `correction::Symbol`: Specifies the type of finite-size correction to apply.
  - `:none`: Returns the raw, uncorrected KBI (`R.kb`). This integral will typically 
    exhibit a drifting tail at long distances due to box size limitations.
  - `:first_order` (default): Applies a geometric, shape-independent first-order 
    dampening window. It maps the arbitrary protein geometry to an effective length 
    scale ``L = 6V/A`` at the final cutoff. A Bartlett-style triangular window 
    ``W(r) = 1 - \\frac{3r}{2L}`` is then applied to the integrand. This perfectly 
    subtracts the leading-order thermodynamic boundary error while dampening 
    long-range statistical noise.

# References
- Krüger, P., & Vlugt, T. J. H. (2018). "Size and shape dependence of finite-volume 
  Kirkwood-Buff integrals." *Physical Review E*, 97(5), 051301. 
  [DOI: 10.1103/PhysRevE.97.051301]

"""
function kbi(R::Result; correction=:first_order)
    dr = R.files[1].options.binstep
    u = units.Angs3tocm3permol
    kb = if correction == :none
        copy(r.kb)
    elseif correction == :first_order
        # Compute L = 6V/A at the final cutoff 
        V_R = sum(R.volume.shell) / R.density.solvent_bulk
        A_R = R.volume.shell[end] / (R.density.solvent_bulk * dr)
        L = 6 * V_R / A_R

        # Leading-order shape-independent correction from Eq. 12 of 
        # Kruger & Vlugt (10.1103/PhysRevE.97.051301)
        W = @. 1 - (3/2) * R.d / L
        u * (1/R.density.solvent_bulk) * cumsum(
            (R.md_count[i] - R.md_count_random[i]) * W[i]
            for i in eachindex(R.d)
        )
    end
    return kb
end

@testitem "kbi" begin

    # needs tests

end