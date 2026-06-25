import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from itertools import groupby
import sys

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_D_wb97x_d_4000_ts.xyz" 
ENERGY_FILE = "CCSD_Combined_Results.txt"
AU_TO_KCAL = 627.509
TS_TO_FS = 6.04721e-16 * 1e15 # Conversion factor

def calculate_dwell_times(df_md, ts_to_fs):
    """Identifies and prints the duration spent in each stabilizing basin."""
    
    # 1. Assign a state to every timestep based on your logic
    # Using 'N' for Neutral (between -6 and 6)
    states = []
    traj_om = df_md['diff_om'].values
    traj_pm = df_md['diff_pm'].values
    
    for om, pm in zip(traj_om, traj_pm):
        if om <= -6 and om < pm:
            states.append("Ortho")
        elif pm <= -6 and pm < om:
            states.append("Para")
        elif om >= 6 and pm >= 6:
            states.append("Meta")
        else:
            states.append("Neutral")

    # 2. Group contiguous states and calculate duration
    # groupby returns (key, group_iterator)
    dwell_data = []
    for state, group in groupby(states):
        count = len(list(group))
        duration_fs = count * ts_to_fs
        dwell_data.append((state, duration_fs))

    # 3. Print the results chronologically
    print(f"\n{'='*40}")
    print(f"{'CHRONOLOGICAL DWELL TIMES':^40}")
    print(f"{'='*40}")
    
    counts = {"Ortho": 0, "Para": 0, "Meta": 0, "Neutral": 0}
    ordinal = ["first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth"]
    
    for state, duration in dwell_data:
        if state == "Neutral": continue # Skip printing the transitions if desired
        
        idx = counts[state]
        prefix = ordinal[idx] if idx < len(ordinal) else f"{idx+1}th"
        
        print(f"{prefix:>10} {state:<6}: {duration:>8.2f} fs")
        counts[state] += 1
    
    print(f"{'='*40}\n")



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

def plot_deltaE_timeseries(df_md, file_name="deltaE_vs_time_direction_D.png"):
    """Plot with comparative shading to show the thermodynamic winner."""
    # 1. Styling
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "font.size": 16,
        "mathtext.fontset": "stix",
    })

    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Time conversion
    ts_to_fs = 6.04721e-16 * 1e15 
    
    # Data Windowing
    slice_idx = slice(10, None)
    t_fs = df_md['step'][slice_idx] * ts_to_fs
    traj_om = df_md['diff_om'][slice_idx].values
    traj_pm = df_md['diff_pm'][slice_idx].values
    traj_rel_e = (df_md['e_md'] - df_md['e_md'].min())[slice_idx].values * AU_TO_KCAL
    
    # Colors
    c_ortho = "#0072B2"   # Deep Teal/Blue
    c_para = "#E69F00"    # Forest/Teal
    c_nitro = "#000000"   # Burnt Orange
    c_destab = "#CC79A7"  # Maroon
    
    # 2. Comparative Shading Logic
    
    # ORTHO WINNER: Ortho < -6 AND Ortho < Para
    ortho_wins = (traj_om <= -6) & (traj_om < traj_pm)
    ax.fill_between(t_fs, traj_om, -6, where=ortho_wins, 
                    color=c_ortho, alpha=0.2, interpolate=True, label='Ortho Preferred')

    # PARA WINNER: Para < -6 AND Para < Ortho
    para_wins = (traj_pm <= -6) & (traj_pm < traj_om)
    ax.fill_between(t_fs, traj_pm, -6, where=para_wins, 
                    color=c_para, alpha=0.2, interpolate=True, label='Para Preferred')

    # BOTH DESTABILIZED: Both > 6
    both_bad = (traj_om >= 6) & (traj_pm >= 6)
    ax.fill_between(t_fs, 6, np.maximum(traj_om, traj_pm), where=both_bad, 
                    color=c_destab, alpha=0.15, interpolate=True, label='Meta Preferred')

    # 3. Trajectory Lines
    ax.plot(t_fs, traj_rel_e, color=c_nitro, lw=1.5, alpha=0.7, label="Nitrobenzene Rel. $E$")
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
    
    # Simplified legend to avoid clutter from fill_between labels
    handles, labels = ax.get_legend_handles_labels()
    # Unique labels only
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=16, frameon=True)

    plt.tight_layout()
    plt.savefig(file_name, dpi=600, bbox_inches="tight")
    print(f"Plot saved: {file_name}")
    plt.show()

def main():
    # 1. Load and Clean
    df_md = parse_md_data(MD_FILE)
    df_grid = pd.read_csv(ENERGY_FILE, sep='\s+', skiprows=[1])
    
    for col in ['theta', 'phi', 'Para_E', 'Ortho_E', 'Meta_E']:
        df_grid[col] = pd.to_numeric(df_grid[col], errors='coerce')

    # 2. Interpolate
    df_md = get_interpolated_energies(df_grid, df_md)

    # 3. Calculate and Print Dwell Times
    calculate_dwell_times(df_md, TS_TO_FS)

    # 4. Generate Plot
    plot_deltaE_timeseries(df_md)

if __name__ == "__main__":
    main()