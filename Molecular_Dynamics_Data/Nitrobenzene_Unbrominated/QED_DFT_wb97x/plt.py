import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import sys

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_A_wb97x_4000_ts.xyz" #"md_trajectory.txt"
ENERGY_FILE = "CCSD_Combined_Results.txt"
AU_TO_KCAL = 627.509

def parse_md_data(filename):
    """Parses MD block format to extract Step, Energy, Phi, and Theta."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if "Step" in line:
                parts = line.split()
                # Parse: "Step 0 E=-436... phi=36.879 theta=73.569"
                data.append({
                    'step': int(parts[1]),
                    'e_md': float(parts[2].split('=')[1]),
                    'phi': float(parts[3].split('=')[1]),
                    'theta': float(parts[4].split('=')[1])
                })
    return pd.DataFrame(data)

def plot_deltaE_timeseries(df_md, file_name="deltaE_vs_time.png"):
    """Publication-ready plot with threshold shading."""
    # 1. Styling
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "font.size": 16,
        "mathtext.fontset": "stix",
    })

    fig, ax = plt.subplots(figsize=(7.5, 5))
    
    # Data extraction
    t_values = df_md['step']
    traj1 = df_md['diff_om_grid'] # Ortho-Meta
    traj2 = (df_md['e_md'] - df_md['e_md'].min()) * AU_TO_KCAL # Rel Energy
    
    color1, color2 = "#005F73", "#BB3E03"

    # 2. Shading for stabilization/destabilization
    ax.fill_between(t_values, traj1, -6, where=(traj1 <= -6), color=color1, alpha=0.15, interpolate=True)
    ax.fill_between(t_values, traj1, 6, where=(traj1 >= 6), color=color1, alpha=0.15, interpolate=True)

    # 3. Plots
    ax.plot(t_values, traj1, color=color1, linewidth=2.0, alpha=0.9, label="Ortho vs Meta ($\Delta E$)")
    ax.plot(t_values, traj2, color=color2, linewidth=2.0, alpha=0.9, label="Nitrobenzene Rel. $E$")

    # 4. Ref lines & Formatting
    ax.axhline(0.0, color='black', linestyle='-', linewidth=1.0, alpha=0.6, zorder=0)
    ax.set_xlabel("Time (steps)", fontsize=16, labelpad=8)
    ax.set_ylabel(r"$\Delta E$ ($\mathrm{kcal \cdot mol^{-1}}$)", fontsize=16, labelpad=8)
    ax.set_ylim(-25, 45)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle=':', alpha=0.3, color='gray')
    ax.legend(loc='best', fontsize=14, frameon=True, edgecolor='0.8')

    plt.tight_layout()
    plt.savefig(file_name, dpi=600, bbox_inches="tight")
    print(f"Plot saved: {file_name}")
    plt.show()

def main():
    # 1. Load Data
    try:
        df_md = parse_md_data(MD_FILE)
        df_grid = pd.read_csv(ENERGY_FILE, delim_whitespace=True, skiprows=[1])
    except Exception as e:
        print(f"Error loading files: {e}")
        return

    # 2. Match MD trajectory to nearest grid point
    # Uses KDTree for efficiency
    tree = KDTree(df_grid[['theta', 'phi']].values)
    _, idx = tree.query(df_md[['theta', 'phi']].values)
    
    df_md['ortho_e_grid'] = df_grid.loc[idx, 'Ortho_E'].values
    df_md['meta_e_grid'] = df_grid.loc[idx, 'Meta_E'].values
    df_md['diff_om_grid'] = (df_md['ortho_e_grid'] - df_md['meta_e_grid']) * AU_TO_KCAL

    # 3. Generate Plot
    plot_deltaE_timeseries(df_md)

if __name__ == "__main__":
    main()
