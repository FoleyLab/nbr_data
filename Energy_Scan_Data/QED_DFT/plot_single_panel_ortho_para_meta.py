import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---
FILE_NAME = "CCSD_Combined_Results.txt"
DPI = 300  
AU_TO_KCAL = 627.509

def create_single_map(P, T, data_grid, title, filename, label, vmin=None, vmax=None):
    """Helper to generate a single publication-ready plot."""
    fig, ax = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    
    # Use a divergent colormap centered at 0
    # Determine symmetric limits for the colorbar
    if vmax==None and vmin==None:
       vmax = np.abs(data_grid).max()
    
       im = ax.pcolormesh(P, T, data_grid, shading='gouraud', cmap='RdBu_r', 
                           vmin=-vmax, vmax=vmax)

    else:
       im = ax.pcolormesh(P, T, data_grid, shading='gouraud', cmap='RdBu_r',
                           vmin=vmin, vmax=vmax)
    
    ax.set_title(title, fontsize=16, fontweight='bold', pad=15)
    ax.set_xlabel(r'Azimuthal Angle $\phi$ ($^{\circ}$)', fontsize=14)
    ax.set_ylabel(r'Polar Angle $\theta$ ($^{\circ}$)', fontsize=14)
    
    cbar = fig.colorbar(im, ax=ax, label=f'Energy Difference ({label})')
    cbar.ax.tick_params(labelsize=12)
    
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 45, 90, 135, 180])
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.savefig(filename, dpi=DPI)
    print(f"Saved: {filename}")
    plt.close(fig)

def plot_energy_diffs():
    # 1. Load the data
    data = np.genfromtxt(FILE_NAME, skip_header=2)
    theta, phi = data[:, 0], data[:, 1]
    e_para, e_ortho, e_meta = data[:, 5], data[:, 6], data[:, 7]

    # 2. Reshape
    num_t = len(np.unique(theta))
    num_p = len(np.unique(phi))
    T = theta.reshape(num_t, num_p)
    P = phi.reshape(num_t, num_p)

    # 3. Calculate Differences
    diff_om = (e_ortho - e_meta).reshape(num_t, num_p) * AU_TO_KCAL
    diff_pm = (e_para - e_meta).reshape(num_t, num_p) * AU_TO_KCAL

    # 4. Generate Separate Plots
    create_single_map(P, T, diff_om, r'$\Delta E$ (Ortho $-$ Meta)', 
                      "ortho_meta_diff.png", "kcal/mol")
    
    create_single_map(P, T, diff_pm, r'$\Delta E$ (Para $-$ Meta)', 
                      "para_meta_diff.png", "kcal/mol") #, -15, 27)

if __name__ == "__main__":
    plot_energy_diffs()
