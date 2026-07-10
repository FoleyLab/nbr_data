import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# --- Configuration ---
FILE_NAME ="isomer_energy_decomposition.dat" #"CCSD_Combined_Results.txt"
DPI = 300  # High resolution for publication
AU_TO_KCAL = 627.509

def plot_energy_diffs(component="el"):
    # 1. Load the data
    # Col 0: Theta, Col 1: Phi, Col 5: Para, Col 6: Ortho, Col 7: Meta
    data = np.genfromtxt(FILE_NAME, skip_header=2)

    # column headers
    # theta       phi           Ex           Ey           Ez            Para_E_el            Para_E_ph           Para_E_blc           Para_E_dse           Ortho_E_el           Ortho_E_ph          Ortho_E_blc          Ortho_E_dse            Meta_E_el            Meta_E_ph           Meta_E_blc           Meta_E_dse
    column_dict = {
        "theta": 0,
        "phi": 1,
        "para_E_el": 5,
        "para_E_ph": 6,
        "para_E_blc": 7,
        "para_E_dse": 8,
        "ortho_E_el": 9,
        "ortho_E_ph": 10,
        "ortho_E_blc": 11,
        "ortho_E_dse": 12,
        "meta_E_el": 13,
        "meta_E_ph": 14,
        "meta_E_blc": 15,
        "meta_E_dse": 16
    }

    theta = data[:, column_dict["theta"]]
    phi = data[:, column_dict["phi"]]
    if component == "el":
        e_para = data[:, column_dict["para_E_el"]]
        e_ortho = data[:, column_dict["ortho_E_el"]]
        e_meta = data[:, column_dict["meta_E_el"]]
        legend_label = "Electronic Energy Difference (kcal / mol)"
    elif component == "ph":
        e_para = data[:, column_dict["para_E_ph"]]
        e_ortho = data[:, column_dict["ortho_E_ph"]]
        e_meta = data[:, column_dict["meta_E_ph"]]
        legend_label = "Photonic Energy Difference (kcal / mol)"
    elif component == "blc":
        e_para = data[:, column_dict["para_E_blc"]]
        e_ortho = data[:, column_dict["ortho_E_blc"]]
        e_meta = data[:, column_dict["meta_E_blc"]]
        legend_label = "Bilinear Coupling Energy Difference (kcal / mol)"
    elif component == "dse":
        e_para = data[:, column_dict["para_E_dse"]]
        e_ortho = data[:, column_dict["ortho_E_dse"]]
        e_meta = data[:, column_dict["meta_E_dse"]]
        legend_label = "Dipole Self Energy Difference (kcal / mol)"

    elif component == "total":
        e_para = data[:, column_dict["para_E_el"]] + data[:, column_dict["para_E_ph"]] + data[:, column_dict["para_E_blc"]] + data[:, column_dict["para_E_dse"]]
        e_ortho = data[:, column_dict["ortho_E_el"]] + data[:, column_dict["ortho_E_ph"]] + data[:, column_dict["ortho_E_blc"]] + data[:, column_dict["ortho_E_dse"]]
        e_meta = data[:, column_dict["meta_E_el"]] + data[:, column_dict["meta_E_ph"]] + data[:, column_dict["meta_E_blc"]] + data[:, column_dict["meta_E_dse"]]
        legend_label = "Total Energy Difference (kcal / mol)"

    # 2. Reshape into 2D grids
    # We know num_theta_vals = 21, num_phi_vals = 21
    num_t = len(np.unique(theta))
    num_p = len(np.unique(phi))
    
    T = theta.reshape(num_t, num_p)
    P = phi.reshape(num_t, num_p)
    
    # Calculate Differences (Convert Hartree to kcal/mol for better readability if desired)
    # 1 Hartree = 627.509 kcal/mol. Let's stay in Hartree but use scientific notation.
    diff_ortho_meta = (e_ortho - e_meta).reshape(num_t, num_p) * AU_TO_KCAL
    diff_para_meta = (e_para - e_meta).reshape(num_t, num_p) * AU_TO_KCAL

    # 3. Setup Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    
    # Set a divergent colormap (RdBu) centered at 0 or a sequential one
    # Use 'RdBu_r' to show relative stabilization/destabilization
    cmap = 'RdBu_r' 

    # Plot 1: Ortho - Meta
    im1 = ax1.pcolormesh(P, T, diff_ortho_meta, shading='gouraud', cmap=cmap)
    ax1.set_title(r'$\Delta E$ (Ortho $-$ Meta)', fontsize=14)
    ax1.set_xlabel(r'Azimuthal Angle $\phi$ (deg)', fontsize=12)
    ax1.set_ylabel(r'Polar Angle $\theta$ (deg)', fontsize=12)
    fig.colorbar(im1, ax=ax1, label=legend_label)

    # Plot 2: Para - Meta
    im2 = ax2.pcolormesh(P, T, diff_para_meta, shading='gouraud', cmap=cmap)
    ax2.set_title(r'$\Delta E$ (Para $-$ Meta)', fontsize=14)
    ax2.set_xlabel(r'Azimuthal Angle $\phi$ (deg)', fontsize=12)
    ax2.set_ylabel(r'Polar Angle $\theta$ (deg)', fontsize=12)
    fig.colorbar(im2, ax=ax2, label=legend_label)

    # Formatting tweaks
    for ax in [ax1, ax2]:
        ax.set_xticks([0, 90, 180, 270, 360])
        ax.set_yticks([0, 45, 90, 135, 180])
        ax.grid(True, linestyle='--', alpha=0.3)


    plt.suptitle(legend_label, fontsize=16, fontweight='bold')
    
    file_name = f"energy_difference_maps_{component}.png"
    # 4. Save and Show
    plt.savefig(file_name, dpi=DPI)
    print(f"Plot saved as {file_name}")
    plt.show()

if __name__ == "__main__":
    plot_energy_diffs(component="el")  # Change to "ph", "blc", or "dse" for other components
    plot_energy_diffs(component="ph")
    plot_energy_diffs(component="blc")
    plot_energy_diffs(component="dse")
    plot_energy_diffs(component="total")
