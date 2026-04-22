import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---
FILE_NAME = "CCSD_Combined_Results.txt"
DPI = 300
AU_TO_KCAL = 627.509

def plot_relative_energies():
    # 1. Load data
    data = np.genfromtxt(FILE_NAME, skip_header=2)
    theta, phi = data[:, 0], data[:, 1]
    
    # Map columns: Para=5, Ortho=6, Meta=7
    energies = {
        'Ortho': data[:, 6],
        'Meta': data[:, 7],
        'Para': data[:, 5]
    }

    # 2. Reshape and Calculate Relative Energy
    num_t = len(np.unique(theta))
    num_p = len(np.unique(phi))
    T = theta.reshape(num_t, num_p)
    P = phi.reshape(num_t, num_p)

    # 3. Setup Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)
    cmap = 'magma' # Dark to Light (Low to High energy)

    # Global max for consistent color scale across all three
    all_rel_energies = []
    for name, e in energies.items():
        # Subtract reference (theta=0, phi=0 index is 0 in the sorted data)
        ref_e = e[0] 
        rel_e = (e - ref_e) * AU_TO_KCAL
        all_rel_energies.append(rel_e)
    
    global_vmax = np.max(all_rel_energies)

    # 4. Plotting
    for i, (name, e) in enumerate(energies.items()):
        ax = axes[i]
        rel_e = (e - e[0]).reshape(num_t, num_p) * AU_TO_KCAL
        
        im = ax.pcolormesh(P, T, rel_e, shading='gouraud', cmap=cmap, vmin=0, vmax=global_vmax)
        
        ax.set_title(f'{name} Relative Energy', fontsize=18, fontweight='bold', pad=15)
        ax.set_xlabel(r'Azimuthal Angle $\phi$ ($^{\circ}$)', fontsize=16)
        if i == 0:
            ax.set_ylabel(r'Polar Angle $\theta$ ($^{\circ}$)', fontsize=16)
        
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xticks([0, 90, 180, 270, 360])
        ax.set_yticks([0, 45, 90, 135, 180])
        ax.grid(True, linestyle='--', alpha=0.3)

    # Colorbar
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8, label='Relative Energy (kcal/mol)')
    cbar.set_label('Relative Energy (kcal/mol)', fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    plt.savefig("relative_energy_landscapes.png", dpi=DPI)
    print("Saved: relative_energy_landscapes.png")
    plt.show()

if __name__ == "__main__":
    plot_relative_energies()
