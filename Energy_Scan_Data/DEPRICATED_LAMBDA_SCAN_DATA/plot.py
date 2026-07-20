import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1. Load data
# Ensure 'molecular_stabilization_data.csv' is in your working directory
df = pd.read_csv('Direction_A_stabilization_data.csv')

# Constants for coordinate transformations
a0_nm = 0.0529177
C_nm3 = 4 * np.pi * (a0_nm**3)

def l_to_v(l):
    return np.where(l == 0, np.inf, C_nm3 / (l**2))

def v_to_l(v):
    return np.sqrt(C_nm3 / v)

# 2. Plotting setup - Publication Quality
plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif'})
fig, ax1 = plt.subplots(figsize=(10, 7), dpi=300)

# Color-blind/Print friendly colors (Okabe-Ito palette)
color_om = '#0072B2' # Blue
color_pm = '#D55E00' # Vermillion

# 3. Data Plots
ax1.plot(df['lambda_au'], df['DeltaE_OM'], 'o-', color=color_om, 
         linewidth=2.5, markersize=8, label=r'$\Delta E_{OM}$')
ax1.plot(df['lambda_au'], df['DeltaE_PM'], 's-', color=color_pm, 
         linewidth=2.5, markersize=8, label=r'$\Delta E_{PM}$')

# 4. 10 kT stabilization line (-6 kcal/mol)
ax1.axhline(y=-6, color='black', linestyle='--', linewidth=1.5, alpha=0.8)
ax1.text(0.012, -5.5, r'~10 $k_B T$ Stabilization threshold', fontsize=12, fontweight='bold')

# 5. Labels and Styling
ax1.set_xlabel(r'Coupling Strength $\lambda$ (a.u.)', fontsize=16, labelpad=10)
ax1.set_ylabel(r'$\Delta E$ (kcal/mol)', fontsize=16, labelpad=10)
ax1.set_title('Stabilization Energy vs. Cavity Parameters', fontsize=18, pad=25)
ax1.grid(True, which='both', linestyle=':', alpha=0.5)
ax1.set_xlim(0.005, 0.105)

# 6. Secondary Axis (Volume)
# This maps the non-linear volume scale to the top of the plot
secax = ax1.secondary_xaxis('top', functions=(l_to_v, v_to_l))
secax.set_xlabel(r'Mode Volume $V$ (nm$^3$)', fontsize=16, labelpad=15)

# We define specific lambda points to display as Volume ticks for clarity
lambda_tick_locs = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
v_ticks = l_to_v(lambda_tick_locs)
secax.set_ticks(v_ticks)
secax.set_xticklabels([f'{v:.2f}' for v in v_ticks])
#secax.set_ticklabels([f'{v:.1f}' for v in v_ticks])

ax1.legend(loc='upper right', frameon=True, fontsize=14)

plt.tight_layout()
plt.savefig('stabilization_plot.png', bbox_inches='tight')
plt.show()
