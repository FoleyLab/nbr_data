import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# --- Configuration ---
DATA_FILE = "para_meta_refinement_scan.txt" #ortho_meta_refinement_scan.txt"
OUTPUT_IMAGE = "para_meta_difference_surface.png"

# Conversion factor: 1 Hartree = 627.509469 kcal/mol
AU_TO_KCAL_MOL = 627.509469

# --- 1. Load and Parse Data ---
# Read the file line by line to skip the header cleanly
data = []
with open(DATA_FILE, "r") as f:
    for line in f:
        # Skip header lines
        if line.strip().startswith("Theta") or line.strip().startswith("---") or not line.strip():
            continue
        data.append([float(x) for x in line.split()])

data = np.array(data)

# Extract column arrays based on your output format:
# Columns: Theta (0), Phi (1), Ex (2), Ey (3), Ez (4), Ortho_E (5), Meta_E (6)
theta = data[:, 0]
phi = data[:, 1]
ortho_e = data[:, 5]
meta_e = data[:, 6]

# --- 2. Process Energies ---
# Calculate the difference and convert to kcal/mol
energy_diff_au = ortho_e - meta_e
energy_diff_kcal = energy_diff_au * AU_TO_KCAL_MOL

# Relative energy: shift the minimum of the difference to 0 for better contrast
#energy_diff_kcal -= np.nanmin(energy_diff_kcal)

# --- 3. Grid Data for Contour Plotting ---
# Create a dense regular mesh grid from our scattered data coordinates
theta_grid = np.linspace(theta.min(), theta.max(), 200)
phi_grid = np.linspace(phi.min(), phi.max(), 200)
T, P = np.meshgrid(theta_grid, phi_grid)

# Interpolate the kcal/mol differences onto the regular mesh grid
Z = griddata((theta, phi), energy_diff_kcal, (T, P), method='cubic')

# --- 4. Plotting ---
plt.figure(figsize=(8, 6.5), dpi=150)

# Draw the smooth color filled contour map
contour_filled = plt.contourf(T, P, Z, levels=30, cmap="viridis")

# Add a colorbar with labels
cbar = plt.colorbar(contour_filled)
cbar.set_label(r"$\Delta E$ (Para - Meta) [kcal/mol] (Relative)", fontsize=12)

# Overlay subtle line contours for easier structure reading
contours = plt.contour(T, P, Z, levels=10, colors="white", alpha=0.3, linewidths=0.7)
plt.clabel(contours, inline=True, fmt='%1.1f', fontsize=8, colors="white")

# Plot your exact calculated scatter points to check your grid coverage
plt.scatter(theta, phi, color="red", s=3, alpha=0.4, label="Data Points")

# Aesthetics
plt.title(r"Energy Surface Refinement: Para vs Meta", fontsize=14, pad=15)
plt.xlabel(r"$\theta$ (degrees)", fontsize=12)
plt.ylabel(r"$\phi$ (degrees)", fontsize=12)
plt.xlim(theta.min(), theta.max())
plt.ylim(phi.min(), phi.max())
plt.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout()

# Save and Show
plt.savefig(OUTPUT_IMAGE)
plt.show()

print(f"Plot successfully saved to {OUTPUT_IMAGE}")
