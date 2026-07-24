
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load the datasets
df_a = pd.read_csv('Direction_A_stabilization_data.csv')
df_d = pd.read_csv('Direction_D_stabilization_data.csv')

# Constants for coordinate transformations
a0_nm = 0.05291772109
C_nm3 = 4 * np.pi * (a0_nm**3)

def l_to_v(l):
    return np.where(l == 0, np.inf, C_nm3 / (l**2))

def v_to_l(v):
    return np.sqrt(C_nm3 / v)

# Publication Quality Setup
plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif'})
fig, ax1 = plt.subplots(figsize=(11, 7.5), dpi=300)

# Colors
color_a = '#0072B2' # Blue (Color-blind friendly)
color_d = '#D55E00' # Vermillion (Color-blind friendly)

# Plot requested columns: OM from A and PM from D
ax1.plot(df_a['lambda_au'], df_a['DeltaE_OM'], 'o-', color=color_a, 
         linewidth=2.5, markersize=9, 
         label=r'$\Delta E_{ortho-meta}, \; \theta=70^{\circ}, \; \phi=31^{\circ}$')

ax1.plot(df_d['lambda_au'], df_d['DeltaE_PM'], 's-', color=color_d, 
         linewidth=2.5, markersize=9, 
         label=r'$\Delta E_{para-meta}, \; \theta=63^{\circ}, \; \phi=63^{\circ}$')

# 10 kT stabilization line (-6 kcal/mol)
ax1.axhline(y=-6, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax1.text(0.007, -5.3, r'~10 $k_B T$ Stabilization', fontsize=12, fontweight='bold', color='#333333')

# Axes labels and styling
ax1.set_xlabel(r'Coupling Strength $\lambda$ (a.u.)', fontsize=16, labelpad=10)
ax1.set_ylabel(r'$\Delta E$ (kcal/mol)', fontsize=16, labelpad=10)
#ax1.set_title('Comparative Stabilization by Coupling Direction', fontsize=18, pad=30)
ax1.grid(True, which='both', linestyle=':', alpha=0.5)
ax1.set_xlim(0.005, 0.105)

# Secondary Axis (Volume in nm^3)
secax = ax1.secondary_xaxis('top', functions=(l_to_v, v_to_l))
secax.set_xlabel(r'Mode Volume $V$ (nm$^3$)', fontsize=16, labelpad=15)

# Align Volume ticks with linear Lambda ticks
lambda_tick_locs = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
v_ticks = l_to_v(lambda_tick_locs)
secax.set_ticks(v_ticks)

# Use .2f formatting as requested and use set_xticklabels for version compatibility
secax.set_xticklabels([f'{v:.2f}' for v in v_ticks])

ax1.legend(loc='lower left', frameon=True, fontsize=13, shadow=True)

plt.tight_layout()
plt.savefig('combined_stabilization_plot.png', bbox_inches='tight')
plt.show()
