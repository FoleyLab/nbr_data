import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_C_wb97x_d_4000_ts.xyz"
ENERGY_FILE = "CCSD_Combined_Results.txt"
AU_TO_KCAL = 627.509
DPI = 600 # High resolution for publication

def parse_md_data(filename):
    """Parses MD block format and applies theta branch correction."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if "Step" in line:
                parts = line.split()
                raw_theta = float(parts[4].split('=')[1])
                
                # Theta phase correction logic
                corrected_theta = raw_theta
                if corrected_theta > 100:
                    corrected_theta = 180.0 - corrected_theta
                
                data.append({
                    'step': int(parts[1]),
                    'phi': float(parts[3].split('=')[1]),
                    'theta': corrected_theta
                })
    return pd.DataFrame(data)

def generate_overlay_plots():
    # 1. Global Styling
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "font.size": 16,
        "mathtext.fontset": "stix",
    })

    # 2. Load and Prepare Grid Data
    df_grid = pd.read_csv(ENERGY_FILE, sep='\s+', skiprows=[1])
    # Ensure numeric types
    for col in ['theta', 'phi', 'Ortho_E', 'Para_E', 'Meta_E']:
        df_grid[col] = pd.to_numeric(df_grid[col], errors='coerce')
    df_grid = df_grid.dropna()

    # Reshape for pcolormesh
    unique_t = np.unique(df_grid['theta'])
    unique_p = np.unique(df_grid['phi'])
    num_t, num_p = len(unique_t), len(unique_p)
    
    T = df_grid['theta'].values.reshape(num_t, num_p)
    P = df_grid['phi'].values.reshape(num_t, num_p)
    
    # Calculate Difference Surfaces
    diff_om = (df_grid['Ortho_E'] - df_grid['Meta_E']).values.reshape(num_t, num_p) * AU_TO_KCAL
    diff_pm = (df_grid['Para_E'] - df_grid['Meta_E']).values.reshape(num_t, num_p) * AU_TO_KCAL

    # 3. Load MD Trajectory
    df_md = parse_md_data(MD_FILE)

    # 4. Define Plots to Generate
    plots = [
        (diff_om, "Ortho - Meta", "traj_overlay_ortho_meta_direction_C.png"),
        (diff_pm, "Para - Meta", "traj_overlay_para_meta_direction_C.png")
    ]

    for surface, title, filename in plots:
        fig, ax = plt.subplots(figsize=(8, 6.5), constrained_layout=True)
        
        # Determine symmetric color limits for zero-centered RdBu
        vlimit = max(abs(surface.min()), abs(surface.max()))
        
        # Plot PES Background
        im = ax.pcolormesh(P, T, surface, cmap='RdBu_r', shading='gouraud', 
                           vmin=-vlimit, vmax=vlimit)
        
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(r'$\Delta E$ (kcal/mol)', fontsize=16)

        # Plot MD Trajectory
        # Black line with alpha for the path
        ax.plot(df_md['phi'], df_md['theta'], color='black', lw=2.0, alpha=0.5, 
                label='MD Trajectory', zorder=2)
        
        # Start/End points with high-contrast markers
        ax.scatter(df_md['phi'].iloc[0], df_md['theta'].iloc[0], 
                   color='#00FF00', edgecolors='black', s=80, label='Start', zorder=3)
        ax.scatter(df_md['phi'].iloc[-1], df_md['theta'].iloc[-1], 
                   color='#FF00FF', edgecolors='black', s=100, marker='X', label='End', zorder=3)

        # Labels and formatting
        ax.set_title(f'Trajectory on {title} Landscape', fontweight='bold', pad=12)
        ax.set_xlabel(r'Azimuthal Angle $\phi$ (deg.)')
        ax.set_ylabel(r'Polar Angle $\theta$ (deg.)')
        
        ax.set_xticks([0, 90, 180, 270, 360])
        ax.set_yticks([0, 45, 90, 135, 180])
        
        ax.legend(loc='lower right', fontsize=12, framealpha=0.9)
        
        plt.savefig(filename, dpi=DPI)
        print(f"Saved: {filename}")
        plt.show()

if __name__ == "__main__":
    generate_overlay_plots()