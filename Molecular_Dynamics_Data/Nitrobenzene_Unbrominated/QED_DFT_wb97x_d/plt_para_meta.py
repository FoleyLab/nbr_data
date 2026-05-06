import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import sys

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_C_wb97x_d_4000_ts.xyz" 
ENERGY_FILE = "CCSD_Combined_Results.txt"
AU_TO_KCAL = 627.509

def get_interpolated_energies(df_grid, df_md):
    """Interpolates QED-CCSD energies (Ortho, Para, Meta) onto MD trajectory."""
    theta_vals = np.sort(df_grid['theta'].unique())
    phi_vals = np.sort(df_grid['phi'].unique())
    
    # Helper to create interpolator for a specific column
    def make_interp(col_name):
        grid_data = df_grid.pivot(index='theta', columns='phi', values=col_name).values
        return RegularGridInterpolator((theta_vals, phi_vals), grid_data, 
                                       bounds_error=False, fill_value=None)

    interp_ortho = make_interp('Ortho_E')
    interp_para = make_interp('Para_E')
    interp_meta = make_interp('Meta_E')
    
    pts = df_md[['theta', 'phi']].values
    
    # Evaluate and store in MD dataframe
    df_md['ortho_e_interp'] = interp_ortho(pts)
    df_md['para_e_interp'] = interp_para(pts)
    df_md['meta_e_interp'] = interp_meta(pts)
    
    # Calculate Differences
    df_md['diff_om'] = (df_md['ortho_e_interp'] - df_md['meta_e_interp']) * AU_TO_KCAL
    df_md['diff_pm'] = (df_md['para_e_interp'] - df_md['meta_e_interp']) * AU_TO_KCAL
    
    return df_md

def parse_md_data(filename):
    """Parses MD blocks and corrects theta branch/phase ambiguity."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if "Step" in line:
                parts = line.split()
                raw_theta = float(parts[4].split('=')[1])
                
                # Apply 180-theta correction logic
                corrected_theta = raw_theta
                if corrected_theta > 100:
                    corrected_theta = 180.0 - corrected_theta
                
                data.append({
                    'step': int(parts[1]),
                    'e_md': float(parts[2].split('=')[1]),
                    'phi': float(parts[3].split('=')[1]),
                    'theta': corrected_theta
                })
    return pd.DataFrame(data)

def plot_deltaE_timeseries(df_md, file_name="deltaE_vs_time_direction_C.png"):
    """Plot with multiple differences and specific conditional shading."""
    # 1. Styling
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "font.size": 16,
        "mathtext.fontset": "stix",
    })

    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Time conversion: 1 au = 0.02418884 fs (6.047...e-16 is related to bohr/vel, using standard au time)
    ts_to_fs = 6.04721e-16 * 1e15  #0.02418884 
    
    # Data Windowing
    slice_idx = slice(10, None)
    t_fs = df_md['step'][slice_idx] * ts_to_fs
    traj_om = df_md['diff_om'][slice_idx]
    traj_pm = df_md['diff_pm'][slice_idx]
    traj_rel_e = (df_md['e_md'] - df_md['e_md'].min())[slice_idx] * AU_TO_KCAL
    
    # Colors
    c_ortho = "#005F73"  # Deep Teal
    c_para = "#0A9396"   # Forest/Caribbean Green
    c_nitro = "#BB3E03"  # Burnt Orange
    c_destab = "#800000" # Maroon for both > 6

    # 2. Shading Logic
    # Ortho-Meta Shading
    ax.fill_between(t_fs, traj_om, -6, where=(traj_om <= -6), color=c_ortho, alpha=0.15, interpolate=True)
    ax.fill_between(t_fs, traj_om, 6, where=(traj_om >= 6), color=c_destab, alpha=0.15, interpolate=True)
    
    # Para-Meta Shading
    ax.fill_between(t_fs, traj_pm, -6, where=(traj_pm <= -6), color=c_para, alpha=0.15, interpolate=True)
    ax.fill_between(t_fs, traj_pm, 6, where=(traj_pm >= 6), color=c_destab, alpha=0.25, interpolate=True)

    # 3. Trajectory Lines
    ax.plot(t_fs, traj_rel_e, color=c_nitro, lw=1.8, alpha=0.8, label="Nitrobenzene Rel. $E$")
    ax.plot(t_fs, traj_om, color=c_ortho, lw=2.2, label=r"$\Delta E_{ortho-meta}$")
    ax.plot(t_fs, traj_pm, color=c_para, lw=2.2, linestyle='--', label=r"$\Delta E_{para-meta}$")

    # 4. Ref lines & Formatting
    ax.axhline(0.0, color='black', linestyle='-', linewidth=1.2, alpha=0.7, zorder=0)
    ax.axhline(6.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.axhline(-6.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    
    ax.set_xlabel("Time (fs)", labelpad=8)
    ax.set_ylabel(r"$\Delta E$ ($\mathrm{kcal \cdot mol^{-1}}$)", labelpad=8)
    ax.set_ylim(-25, 45)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle=':', alpha=0.2, color='gray')
    ax.legend(loc='upper right', fontsize=12, frameon=True, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(file_name, dpi=600, bbox_inches="tight")
    print(f"Plot saved: {file_name}")
    plt.show()

def main():
    # 1. Load and Clean
    df_md = parse_md_data(MD_FILE)
    df_grid = pd.read_csv(ENERGY_FILE, sep='\s+', skiprows=[1])
    
    # Force numeric for grid
    for col in ['theta', 'phi', 'Para_E', 'Ortho_E', 'Meta_E']:
        df_grid[col] = pd.to_numeric(df_grid[col], errors='coerce')

    # 2. Interpolate Energies
    df_md = get_interpolated_energies(df_grid, df_md)

    # 3. Generate Plot
    plot_deltaE_timeseries(df_md)

if __name__ == "__main__":
    main()