"""
Time series plot: Ortho-Meta AND Para-Meta energy differences + relative
nitrobenzene energy, along a single AIMD trajectory.

Purpose
-------
Despite the "_para_meta" filename, this script plots BOTH
DeltaE_ortho-meta(t) and DeltaE_para-meta(t) (see flagged issues re:
misleading name, same issue as plot_ortho_para_trajectory_overlays.py), interpolated onto
one AIMD trajectory, together with the trajectory's relative electronic
energy. Shaded regions mark where each intermediate is stabilized
(<= -5 kcal/mol) or both intermediates are destabilized (>= +5 kcal/mol)
relative to Meta.

Inputs
------
- MD_FILE: AIMD trajectory, ".xyz + Step header" format (hardcoded to
  direction C here).
- ENERGY_FILE: whitespace-delimited surface grid, columns
      theta  phi  Ex  Ey  Ez  Para_E  Ortho_E  Meta_E
  (Hartree). RESOLVED: was set to a nonexistent "CCSD_Combined_Results.txt";
  confirmed with Jay 2026-07-10 and set to
  "isomer_Nel_49_Nph_10_total_energies.dat" (90x90 grid, Nel=49
  electronic states / Nph=10 photon Fock states). See README's
  "Swapping in new intermediate-energy data" section to point this at a
  different level of theory later.

Output
------
- deltaE_vs_time_direction_C.png

Units
-----
- Angles in degrees; energies converted Hartree -> kcal/mol (AU_TO_KCAL)
  for plotting.
- Time axis: step index -> femtoseconds via ts_to_fs in
  plot_deltaE_timeseries. CONFIRMED with Jay (2026-07-10): each recorded
  MD step is 25 atomic units of time apart, so
  ts_to_fs = 25 * 0.0241888 fs/au = 0.60472 fs/step.

Assumptions / non-obvious logic
--------------------------------
- Theta branch correction (theta -> 180 - theta when raw theta > 100
  degrees): CONFIRMED with Jay (2026-07-10) as an intentional, empirical
  fix for arccos branch jumps in the upstream angle calculation (theta
  is computed there as np.degrees(np.arccos(cos_theta)), whose principal
  range is [0, 180] but which can jump discontinuously between frames
  due to reference-axis sign ambiguity or numerical sensitivity near
  cos_theta = +/-1). Not a statement about surface symmetry. Left
  unchanged; see plot_para_meta_trajectory_overlay.py's module docstring for the full explanation.
- Bilinear interpolation (RegularGridInterpolator, fill_value=None) is
  used to map trajectory (theta, phi) points onto the Ortho/Para/Meta
  surfaces; out-of-grid-bounds points are silently extrapolated rather
  than erroring - see plot_para_meta_timeseries.py's docstring for the same caveat.
- Shading threshold of +/-5 kcal/mol is hardcoded in three places in
  plot_deltaE_timeseries, matching the project's stated classification
  scheme. (Originally hardcoded as +/-6 kcal/mol; corrected to +/-5 per
  Jay's confirmation 2026-07-10.)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import sys

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_A_wb97x_d_4000_ts.xyz"
ENERGY_FILE = "isomer_Nel_49_Nph_10_total_energies.dat"  # fixed 2026-07-10: was a nonexistent placeholder path
AU_TO_KCAL = 627.509  # Hartree -> kcal/mol conversion factor

def get_interpolated_energies(df_grid, df_md):
    """Interpolates QED-CCSD energies (Ortho, Para, Meta) onto MD trajectory.

    Builds one bilinear RegularGridInterpolator per intermediate (Ortho,
    Para, Meta) over the (theta, phi) surface grid, evaluates each at
    every MD trajectory (theta, phi) point, and stores both the raw
    interpolated energies and the two energy differences of interest
    (Ortho-Meta, Para-Meta), converted to kcal/mol.
    """
    theta_vals = np.sort(df_grid['theta'].unique())
    phi_vals = np.sort(df_grid['phi'].unique())

    # Helper to create interpolator for a specific column
    # NOTE: fill_value=None -> extrapolates rather than erroring/NaN-ing
    # for MD points outside the sampled (theta, phi) grid range.
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

    # Calculate Differences (Hartree -> kcal/mol)
    df_md['diff_om'] = (df_md['ortho_e_interp'] - df_md['meta_e_interp']) * AU_TO_KCAL
    df_md['diff_pm'] = (df_md['para_e_interp'] - df_md['meta_e_interp']) * AU_TO_KCAL

    return df_md

def parse_md_data(filename):
    """Parses MD blocks and corrects theta branch/phase ambiguity.

    Extracts (step, energy, phi, theta) from each "Step ..." header
    line, folding theta -> 180 - theta whenever the raw value exceeds
    100 degrees (confirmed, intentional arccos branch-jump correction;
    see module docstring).
    """
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

def plot_deltaE_timeseries(df_md, file_name="deltaE_vs_time_direction_A.png"):
    """Plot with multiple differences and specific conditional shading.

    Plots DeltaE_ortho-meta(t), DeltaE_para-meta(t), and the trajectory's
    relative electronic energy, shading each DeltaE curve's stabilized
    (<= -5 kcal/mol) and destabilized (>= +5 kcal/mol) regions. The first
    10 recorded frames are dropped (equilibration/startup transient).
    """
    # 1. Styling
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "font.size": 16,
        "mathtext.fontset": "stix",
    })

    fig, ax = plt.subplots(figsize=(9, 6))

    # Time conversion: 1 au = 0.02418884 fs; each recorded MD step is 25 au
    # apart (CONFIRMED with Jay, 2026-07-10), so ts_to_fs = 25 * 0.0241888.
    ts_to_fs = 6.04721e-16 * 1e15  # 25 au/step -> fs/step

    # Data Windowing (skip first 10 frames - startup/equilibration window)
    slice_idx = slice(10, None)
    t_fs = df_md['step'][slice_idx] * ts_to_fs
    traj_om = df_md['diff_om'][slice_idx]
    traj_pm = df_md['diff_pm'][slice_idx]
    traj_rel_e = (df_md['e_md'] - df_md['e_md'].min())[slice_idx] * AU_TO_KCAL

    # Colors
    c_ortho = "#005F73"  # Deep Teal
    c_para = "#0A9396"   # Forest/Caribbean Green
    c_nitro = "#BB3E03"  # Burnt Orange
    c_destab = "#800000" # Maroon for both > 5

    # 2. Shading Logic
    # +/-5 kcal/mol classification threshold (corrected 2026-07-10 from an
    # originally-hardcoded +/-6 to match the project's stated definition).
    # Ortho-Meta Shading
    bound = 0
    ax.fill_between(t_fs, traj_om, -bound, where=(traj_om <= -bound), color=c_ortho, alpha=0.15, interpolate=True)
    ax.fill_between(t_fs, traj_om, bound, where=(traj_om >= bound), color=c_destab, alpha=0.15, interpolate=True)

    # Para-Meta Shading
    ax.fill_between(t_fs, traj_pm, -bound, where=(traj_pm <= -bound), color=c_para, alpha=0.15, interpolate=True)
    ax.fill_between(t_fs, traj_pm, bound, where=(traj_pm >= bound), color=c_destab, alpha=0.25, interpolate=True)

    # 3. Trajectory Lines
    ax.plot(t_fs, traj_rel_e, color=c_nitro, lw=1.8, alpha=0.8, label="Nitrobenzene Rel. $E$")
    ax.plot(t_fs, traj_om, color=c_ortho, lw=2.2, label=r"$\Delta E_{ortho-meta}$")
    ax.plot(t_fs, traj_pm, color=c_para, lw=2.2, linestyle='--', label=r"$\Delta E_{para-meta}$")

    # 4. Ref lines & Formatting
    ax.axhline(0.0, color='black', linestyle='-', linewidth=1.2, alpha=0.7, zorder=0)
    #ax.axhline(5.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    #ax.axhline(-5.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)

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
    """Load MD trajectory + energy grid, interpolate, and plot the time series."""
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
