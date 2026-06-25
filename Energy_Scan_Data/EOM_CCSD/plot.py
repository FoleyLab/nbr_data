import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load data from CSV
df = pd.read_csv("polariton_energy_scan.csv")

# 2. Create an isolated column for the Total Energy Shift
# Subtract the minimum total energy so the variation is clear relative to baseline
df['E_total_shift'] = df['E_total'] - df['E_el_uncoupled']	

E_ccsd = df['E_el_uncoupled'].min()
max_e_var = df['E_total_shift'].max()

# 3. Define the list of columns to plot, their titles, and preferred colormaps
plot_configs = [
    {"col": "E_el",         "title": "Electronic Energy Contribution", "cmap": "rocket"},
    {"col": "E_ph",         "title": "Photon Energy Contribution",     "cmap": "mako"},
    {"col": "E_blc",        "title": "Bilinear Light-Matter Coupling", "cmap": "coolwarm"}, # Diverging cmap for pos/neg
    {"col": "E_dse",        "title": "Dipole Self-Energy (DSE)",       "cmap": "crest"},
    {"col": "E_total_shift","title": f"Total Energy Variation (Base: {E_ccsd:.6f} Ha, Max dE: {max_e_var:6f})", "cmap": "inferno"}
]

# 4. Generate the plots
for config in plot_configs:
    plt.figure(figsize=(8, 6))
    
    # Pivot data from long format to a 2D Grid
    grid_df = df.pivot(index="theta", columns="phi", values=config["col"])
    
    # Sort index so theta runs logically from 0 to 180 top-to-bottom or bottom-to-top
    grid_df = grid_df.sort_index(ascending=False) 
    
    # Create the heatmap
    # robust=True automatically trims extreme outliers from color scaling if needed
    ax = sns.heatmap(
        grid_df, 
        cmap=config["cmap"], 
        cbar_kws={'label': 'Energy (Hartree)'},
        robust=True 
    )
    
    # Clean up axis labels so they don't crowd the plot
    # Show labels every 5th tick mark
    x_ticks = [int(float(label.get_text())) for label in ax.get_xticklabels()]
    y_ticks = [int(float(label.get_text())) for label in ax.get_yticklabels()]
    ax.set_xticks(ax.get_xticks()[::5])
    ax.set_xticklabels(x_ticks[::5])
    ax.set_yticks(ax.get_yticks()[::5])
    ax.set_yticklabels(y_ticks[::5])
    
    plt.title(config["title"], fontsize=14, pad=15)
    plt.xlabel(r"$\phi$ (Degrees)", fontsize=12)
    plt.ylabel(r"$\theta$ (Degrees)", fontsize=12)
    
    plt.tight_layout()
    
    # Save each figure separately
    filename = f"heatmap_{config['col']}.png"
    plt.savefig(filename, dpi=300)
    print(f"Generated plot: {filename}")
    
    plt.show()
