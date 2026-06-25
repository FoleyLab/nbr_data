import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---
FILE_NAME = "intermediate_scans.csv"
DPI = 350  
AU_TO_KCAL = 627.509


def create_single_map(P, T, data_grid, title, filename, label, vmin=None, vmax=None):
    """Helper to generate a single publication-ready plot with larger fonts."""
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
    fs_title = 20
    fs_labels = 18
    fs_ticks = 16
    
    if vmin is None and vmax is None:
        vmax = np.abs(data_grid).max()
        im = ax.pcolormesh(P, T, data_grid, shading='gouraud', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    else:
        im = ax.pcolormesh(P, T, data_grid, shading='gouraud', cmap='RdBu_r', vmin=vmin, vmax=vmax)

    ax.set_title(title, fontsize=fs_title, fontweight='bold', pad=20)
    ax.set_xlabel(r'Azimuthal Angle $\phi$ (deg.)', fontsize=fs_labels)
    ax.set_ylabel(r'Polar Angle $\theta$ (deg.)', fontsize=fs_labels)
    
    ax.tick_params(axis='both', which='major', labelsize=fs_ticks)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 45, 90, 135, 180])
    
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'Energy Difference ({label})', fontsize=fs_labels)
    cbar.ax.tick_params(labelsize=fs_ticks)
    
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.savefig(filename, dpi=DPI)
    print(f"Saved: {filename}")
    plt.close(fig)


def plot_energy_diffs():
    # 1. Load the existing data
    data = np.genfromtxt(FILE_NAME, skip_header=1, delimiter=",")
    theta_raw, phi_raw = data[:, 0], data[:, 1]
    e_ortho_raw, e_meta_raw, e_para_raw = data[:, 3], data[:, 5], data[:, 7]
    
    # Map the raw data into a dictionary for fast coordinate lookup
    # We round to 2 decimal places to prevent floating point mismatch issues
    lookup = {}
    for t, p, eo, em, ep in zip(theta_raw, phi_raw, e_ortho_raw, e_meta_raw, e_para_raw):
        lookup[(round(t, 2), round(p, 2))] = (eo, em, ep)

    # 2. Define the FULL grid boundaries you want to plot
    # Automatically extracts unique step sizes used in your data
    d_theta = np.diff(np.unique(theta_raw)).min()
    d_phi = np.diff(np.unique(phi_raw)).min()
    
    full_theta = np.arange(0, 180 + d_theta/2, d_theta)
    full_phi = np.arange(0, 360 + d_phi/2, d_phi)
    
    num_t = len(full_theta)
    num_p = len(full_phi)
    
    T, P = np.meshgrid(full_theta, full_phi, indexing='ij')
    
    grid_ortho = np.zeros((num_t, num_p))
    grid_meta = np.zeros((num_t, num_p))
    grid_para = np.zeros((num_t, num_p))

    # 3. Populate full grid utilizing vec -> -vec symmetry where needed
    for t_idx, t in enumerate(full_theta):
        for p_idx, p in enumerate(full_phi):
            t_key, p_key = round(t, 2), round(p, 2)
            
            # Scenario A: The point already exists in the CSV data
            if (t_key, p_key) in lookup:
                grid_ortho[t_idx, p_idx], grid_meta[t_idx, p_idx], grid_para[t_idx, p_idx] = lookup[(t_key, p_key)]
            
            # Scenario B: Point is missing -> Find its inverted partner (-vec)
            else:
                t_inv = round(180.0 - t, 2)
                p_inv = round((p + 180.0) % 360.0, 2)
                
                # Handle edge-case wrapping where 360.0 degrees rolls back to 0.0
                if p_inv == 360.0: 
                    p_inv = 0.0
                
                if (t_inv, p_inv) in lookup:
                    grid_ortho[t_idx, p_idx], grid_meta[t_idx, p_idx], grid_para[t_idx, p_idx] = lookup[(t_inv, p_inv)]
                else:
                    # Fallback if neither point nor its inverse exists (should not happen if data is 1/2 complete)
                    grid_ortho[t_idx, p_idx] = np.nan
                    grid_meta[t_idx, p_idx] = np.nan
                    grid_para[t_idx, p_idx] = np.nan

    # 4. Calculate Differences in kcal/mol
    diff_om = (grid_ortho - grid_meta) * AU_TO_KCAL
    diff_pm = (grid_para - grid_meta) * AU_TO_KCAL

    # 5. Generate Separate Plots
    create_single_map(P, T, diff_om, r'$\Delta E$ (Ortho $-$ Meta)', 
                      "ortho_meta_diff_QED_CCSD_22.png", "kcal/mol", vmin=-10, vmax=10)
    
    create_single_map(P, T, diff_pm, r'$\Delta E$ (Para $-$ Meta)', 
                      "para_meta_diff_QED_CCSD_22.png", "kcal/mol", vmin=-10, vmax=10)


if __name__ == "__main__":
    plot_energy_diffs()
