import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import sys

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_B_wb97x_4000_ts.xyz" #"md_trajectory.txt"
ENERGY_FILE = "CCSD_Combined_Results.txt"
AU_TO_KCAL = 627.509

from scipy.interpolate import RegularGridInterpolator

def get_interpolated_energies(df_grid, df_md):
    """
    Interpolates QED-CCSD energies onto MD trajectory points.
    """
    # 1. Reshape the coarse grid into regular 2D arrays
    # Assuming your grid is theta (num_t) x phi (num_p)
    theta_vals = np.sort(df_grid['theta'].unique())
    phi_vals = np.sort(df_grid['phi'].unique())
    
    # Create the interpolator function for Ortho and Meta
    # 'linear' is usually best for energy surfaces to avoid artificial oscillations
    interp_ortho = RegularGridInterpolator((theta_vals, phi_vals), 
                                           df_grid.pivot(index='theta', columns='phi', values='Ortho_E').values,
                                           bounds_error=False, fill_value=None)
    
    interp_meta = RegularGridInterpolator((theta_vals, phi_vals), 
                                          df_grid.pivot(index='theta', columns='phi', values='Meta_E').values,
                                          bounds_error=False, fill_value=None)
    
    # 2. Extract MD coordinates as a list of (theta, phi) pairs
    pts = df_md[['theta', 'phi']].values
    
    # 3. Evaluate the interpolator at these points
    df_md['ortho_e_interp'] = interp_ortho(pts)
    df_md['meta_e_interp'] = interp_meta(pts)
    
    # 4. Calculate smoothed difference
    df_md['diff_om_grid'] = (df_md['ortho_e_interp'] - df_md['meta_e_interp']) * AU_TO_KCAL
    
    return df_md


def parse_md_data(filename):
    """Parses MD block format to extract Step, Energy, Phi, and Theta."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if "Step" in line:
                parts = line.split()
                
                # Extract raw theta
                raw_theta = float(parts[4].split('=')[1])
                
                # Apply correction logic
                corrected_theta = raw_theta
                if corrected_theta > 100:
                    print("correcting theta > 100")
                    corrected_theta = 180.0 - corrected_theta
                # Parse: "Step 0 E=-436... phi=36.879 theta=73.569"
                data.append({
                    'step': int(parts[1]),
                    'e_md': float(parts[2].split('=')[1]),
                    'phi': float(parts[3].split('=')[1]),
                    'theta': corrected_theta
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
    # convert between time step number and femtoseconds
    ts_to_fs = 6.04721e-16 * 1e15 # 1 au of time in femtoseconds
    # Data extraction
    t_values = df_md['step'][10:] * ts_to_fs # Example: first 100 steps
    traj1 = df_md['diff_om_grid'][10:] # Ortho-Meta
    traj2 = (df_md['e_md'] - df_md['e_md'].min())[10:] * AU_TO_KCAL # Rel Energy
    
    # maroon hex color: #800000
    color1, color2, color3 = "#005F73", "#BB3E03", "#800000"

    # 2. Shading for stabilization/destabilization
    ax.fill_between(t_values, traj1, -6, where=(traj1 <= -6), color=color1, alpha=0.15, interpolate=True)
    ax.fill_between(t_values, traj1, 6, where=(traj1 >= 6), color=color3, alpha=0.15, interpolate=True)

    # 3. Plots
    #ax.plot(t_values, traj1, color=color1, linewidth=2.0, alpha=0.9, label="Ortho vs Meta ($\Delta E$)")
    ax.plot(t_values, traj2, color=color2, linewidth=2.0, alpha=0.9, label="Nitrobenzene Relative Energy (kcal / mol)")
    ax.plot(t_values, traj1, color=color1, linewidth=2.0, alpha=0.9, label="$\Delta E_{ortho-meta}$ (kcal / mol)")
    # 4. Ref lines & Formatting
    ax.axhline(0.0, color='black', linestyle='-', linewidth=1.0, alpha=0.6, zorder=0)
    ax.set_xlabel("Time (femtoseconds)", fontsize=16, labelpad=8)
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

#def main():
#    # 1. Load Data
#    try:
#        df_md = parse_md_data(MD_FILE)
#        df_grid = pd.read_csv(ENERGY_FILE, sep='\s+', skiprows=[1])
#    except Exception as e:
#        print(f"Error loading files: {e}")
#        return#
#
#    # 2. Match MD trajectory to nearest grid point
#    # Uses KDTree for efficiency
#    tree = KDTree(df_grid[['theta', 'phi']].values)
#    _, idx = tree.query(df_md[['theta', 'phi']].values)
#    
#    df_md['ortho_e_grid'] = df_grid.loc[idx, 'Ortho_E'].values
#    df_md['meta_e_grid'] = df_grid.loc[idx, 'Meta_E'].values
#    df_md['diff_om_grid'] = (df_md['ortho_e_grid'] - df_md['meta_e_grid']) * AU_TO_KCAL#
#
#    # 3. Generate Plot
#    plot_deltaE_timeseries(df_md)
def main():
    # 1. Load and Clean
    df_md = parse_md_data(MD_FILE)
    df_grid = pd.read_csv(ENERGY_FILE, sep='\s+', skiprows=[1])
    
    # Force numeric for grid
    for col in ['theta', 'phi', 'Ortho_E', 'Meta_E']:
        df_grid[col] = pd.to_numeric(df_grid[col], errors='coerce')

    # 2. Interpolate Energies (Replaces the KDTree block)
    df_md = get_interpolated_energies(df_grid, df_md)

    # 3. Generate Plot
    plot_deltaE_timeseries(df_md)

if __name__ == "__main__":
    main()
