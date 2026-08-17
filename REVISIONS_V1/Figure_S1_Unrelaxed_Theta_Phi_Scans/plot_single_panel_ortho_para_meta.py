import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---
DPI = 100  
AU_TO_KCAL = 627.509

SHADING_OP = "gouraud" #< smooth

# Configuration for different datasets
DATASETS = [
    {"file": "unrelaxed_QED_CCSD_theta_phi_scan.txt", "suffix": "QED\_CCSD"},
    {"file": "unrelaxed_pQED_Nel_49_Nph_10_theta_phi_scan.txt", "suffix": "pQED"},
    {"file": "unrelaxed_pQED_Nel_49_Nph_10_CS_theta_phi_scan.txt", "suffix": "pQED\_CS"},
    {"file": "unrelaxed_qeddft_wb97x_d_theta_phi_scan.txt", "suffix": "qeddf"}
]

def create_single_map(P, T, data_grid, title, filename, label, vmin=None, vmax=None):
    """Helper to generate a single publication-ready plot with larger fonts."""
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
    # 1. Define font sizes
    fs_title = 20
    fs_labels = 18
    fs_ticks = 16
    
    # Colormap scaling
    if vmin==None and vmax==None:
        vmax = np.abs(data_grid).max()
        im = ax.pcolormesh(P, T, data_grid, shading=SHADING_OP, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    else:
        im = ax.pcolormesh(P, T, data_grid, shading=SHADING_OP, cmap='RdBu_r',vmin=vmin, vmax=vmax)

    # 2. Title and Axes labels
    ax.set_title(title, fontsize=fs_title, fontweight='bold', pad=20)
    ax.set_xlabel(r'Azimuthal Angle $\phi$ (deg.)', fontsize=fs_labels)
    ax.set_ylabel(r'Polar Angle $\theta$ (deg.)', fontsize=fs_labels)
    
    # 3. Ticks (for both x and y)
    ax.tick_params(axis='both', which='major', labelsize=fs_ticks)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 45, 90, 135, 180])
    
    # 4. Colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'Energy Difference ({label})', fontsize=fs_labels)
    cbar.ax.tick_params(labelsize=fs_ticks)
    
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.savefig(filename, dpi=DPI)
    print(f"Saved: {filename}")
    plt.close(fig)


def plot_energy_diffs():
    for dataset in DATASETS:
        file_name = dataset["file"]
        suffix = dataset["suffix"]
        
        # 1. Load the data
        data = np.genfromtxt(file_name, skip_header=2)
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
                          f"unrelaxed_ortho_meta_{suffix}.png", "kcal/mol", vmin=-8, vmax=8)
        
        create_single_map(P, T, diff_pm, r'$\Delta E$ (Para $-$ Meta)', 
                          f"unrelaxed_para_meta_{suffix}.png", "kcal/mol", vmin=-8, vmax=8)


if __name__ == "__main__":
    plot_energy_diffs()
