# Activate environment in current directory
import Pkg;
Pkg.activate(".");

# Run this once, to install necessary packages:
# Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings"])

# Load packages
using ComplexMixtures
using PDBTools
using Plots, Plots.Measures
using LaTeXStrings
using Statistics

# Load PDB file of the system
atoms = read_pdb("./system.pdb")

# The results for glycerol were computed in script1.jl
glyc_results = load("./glyc50_results.json")

# Compute the MDDF of water relative to the protein, from the same trajectory
protein = select(atoms, "protein")
water = select(atoms, "water")
solute = AtomSelection(protein, nmols=1)
solvent = AtomSelection(water, natomspermol=3)
trajectory_file = "./glyc50_traj.dcd"
water_results = mddf(trajectory_file, solute, solvent, Options(bulk_range=(10.0, 15.0)))
save(water_results, "water_results.json")
println("Results saved to water_results.json")

#
# KB integrals (converted from cm³/mol to L/mol), and the bulk concentration of glycerol (mol/L)
#
d = glyc_results.d
G_pc = glyc_results.kb / 1000 # protein-glycerol
G_pw = water_results.kb / 1000 # protein-water
ρ_c = overview(glyc_results).density.solvent_bulk

# Preferential interaction parameter, as a function of the distance
Γ = ρ_c .* (G_pc .- G_pw)

# Converged values: average in the region where the KB integrals are stable
converged = findall(r -> 10.0 <= r <= 15.0, d)
Γ_conv = mean(Γ[converged])
ΔG_conv = mean(G_pc[converged] .- G_pw[converged])

# Derivative of the chemical potential of the protein with respect to the glycerol
# concentration, assuming that the water-glycerol mixture is ideal (kcal mol⁻¹ M⁻¹)
T = 298.15 # K
RT = 1.987204e-3 * T # kcal/mol
dμ_dρ_ideal = -RT * ΔG_conv

println("Bulk concentration of glycerol: $(round(ρ_c; digits=2)) mol/L")
println("G_pc - G_pw = $(round(ΔG_conv; digits=1)) L/mol")
println("Preferential interaction parameter, Γ = $(round(Γ_conv; digits=1))")
println("∂μ_p/∂ln(a_c) = -RTΓ = $(round(-RT * Γ_conv; digits=1)) kcal/mol")
println("∂μ_p/∂ρ_c (ideal mixture) = $(round(dμ_dρ_ideal; digits=1)) kcal mol⁻¹ M⁻¹")

#
# Plots
#
Plots.default(
    fontfamily="Computer Modern",
    linewidth=2,
    framestyle=:box,
    label=nothing,
    grid=false
)
plot(layout=(1, 2))
plot!(d, G_pc, label="Glycerol", subplot=1)
plot!(d, G_pw, label="Water", subplot=1)
plot!(xlabel=L"r/\AA", ylabel=L"G_{pi}/\mathrm{L~mol^{-1}}", subplot=1)
plot!(d, Γ, subplot=2)
hline!([0], linestyle=:dash, linecolor=:gray, subplot=2)
plot!(xlabel=L"r/\AA", ylabel=L"\Gamma_{pc}", subplot=2)
plot!(size=(800, 300), margin=4mm)
savefig("./preferential_interaction.png")
println("Created plot preferential_interaction.png")
