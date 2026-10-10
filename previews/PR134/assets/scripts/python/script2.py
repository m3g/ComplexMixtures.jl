import ComplexMixtures as cm
import matplotlib.pyplot as plt

# Load the actual results obtained with the complete simulation:
glyc_results = cm.load("./glyc50.json")
water_results = cm.load("./water.json")

# Plot MDDF and KB
fig, axs = plt.subplots(2)
axs[0].plot(glyc_results.d, glyc_results.mddf, label="Glycerol")
axs[0].plot(water_results.d, water_results.mddf, label="Water")
axs[0].set(ylabel="MDDF")

# Plot KB integral
axs[1].plot(glyc_results.d, glyc_results.kb)
axs[1].plot(water_results.d, water_results.kb)
axs[1].set(xlabel="distance / Angs", ylabel="KB integral")
plt.tight_layout()

plt.savefig("mddf_kb.png")