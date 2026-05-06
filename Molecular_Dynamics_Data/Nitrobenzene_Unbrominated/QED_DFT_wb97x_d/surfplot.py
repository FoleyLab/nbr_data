import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
from scipy.spatial import KDTree

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_D_wb97x_d_4000_ts.xyz"  # Your MD output file
ENERGY_FILE = "CCSD_Combined_Results.txt"
AU_TO_KCAL = 627.509

def parse_md_data(filename):
    """Parses your specific MD block format."""
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "Step" in line:
                
                
                parts = line.split()

                # Extract raw theta
                raw_theta = float(parts[4].split('=')[1])
                corrected_theta = raw_theta
                if corrected_theta > 100:
                    print("correcting theta > 100")
                    corrected_theta = 180.0 - corrected_theta

                step = int(parts[1])
                e = float(parts[2].split('=')[1])
                phi = float(parts[3].split('=')[1])
                #theta = float(parts[4].split('=')[1])
                data.append({'step': step, 'e_md': e, 'theta': corrected_theta, 'phi': phi})
    return pd.DataFrame(data)

# 1. Load Data
df_md = parse_md_data(MD_FILE)
df_grid = pd.read_csv(ENERGY_FILE, sep='\s+', skiprows=[1])

# 2. Map MD coordinates to Grid points (Nearest Neighbor)
# Create a KDTree for fast spatial searching
tree = KDTree(df_grid[['theta', 'phi']].values)
_, idx = tree.query(df_md[['theta', 'phi']].values)

# Get the corresponding Para/Meta energies from the grid
df_md['para_e_grid'] = df_grid.loc[idx, 'Para_E'].values
df_md['meta_e_grid'] = df_grid.loc[idx, 'Meta_E'].values
df_md['diff_om_grid'] = (df_md['para_e_grid'] - df_md['meta_e_grid']) * AU_TO_KCAL

# --- PLOTTING ---

# FIGURE 1: Trajectory on Energy Landscape
plt.figure(figsize=(8, 6))
# set global font properties for publication-quality
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 16,
    "mathtext.fontset": "stix",
})
# Assume we have the grid data from before (reshape)
num_t = len(np.unique(df_grid['theta']))
num_p = len(np.unique(df_grid['phi']))
T = df_grid['theta'].values.reshape(num_t, num_p)
P = df_grid['phi'].values.reshape(num_t, num_p)
diff_om = (df_grid['Para_E'] - df_grid['Meta_E']).values.reshape(num_t, num_p) * AU_TO_KCAL

plt.pcolormesh(P, T, diff_om, cmap='RdBu_r', shading='gouraud')
plt.colorbar(label='Para - Meta Energy (kcal/mol)')
plt.plot(df_md['phi'], df_md['theta'], color='black', lw=2.5, alpha=0.6, label='MD Trajectory')
plt.scatter(df_md['phi'].iloc[0], df_md['theta'].iloc[0], color='green', marker='o', label='Start')
plt.scatter(df_md['phi'].iloc[-1], df_md['theta'].iloc[-1], color='red', marker='X', label='End')
plt.xlabel(r'Azimuthal Angle, $\phi$ (deg.)')
plt.ylabel(r'Polar Angle, $\theta$ (deg.)')
plt.title('MD Trajectory over Para-Meta PES')
plt.legend()
plt.savefig("traj_overlay_direction_D.png", dpi=300)


plt.show()
