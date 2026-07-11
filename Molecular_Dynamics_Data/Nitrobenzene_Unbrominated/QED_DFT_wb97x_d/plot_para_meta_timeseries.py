"""
Time series plot: Para-Meta energy difference + relative nitrobenzene
energy, along a single AIMD trajectory.

Purpose
-------
Interpolates the static Para_E/Meta_E surface onto an AIMD trajectory's
(theta, phi) path (bilinear interpolation, not nearest-neighbor - see
below), then plots DeltaE_para-meta(t) alongside the trajectory's own
relative electronic energy (a proxy for energy-conservation quality),
with shaded regions marking where Para is stabilized/destabilized by
more than 6 kcal/mol relative to Meta.

Inputs
------
- MD_FILE: AIMD trajectory in the project's ".xyz + Step header" format.
  RESOLVED (was a bug): originally read
  "nitrobenzene_direction__wb97x_d_4000_ts.xyz" (blank direction letter,
  matched no file). Confirmed with Jay 2026-07-10 and set to direction B
  to match this script's default output filename
  ("deltaE_vs_time_direction_B.png"), the one direction with no existing
  output figure.
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
- deltaE_vs_time_direction_B.png (default; pass a different file_name to
  plot_deltaE_timeseries to change it).

Units
-----
- Angles in degrees; energies in Hartree from the input files, converted
  to kcal/mol (AU_TO_KCAL) for plotting.
- Time axis: trajectory step index converted to femtoseconds via
  ts_to_fs inside plot_deltaE_timeseries. CONFIRMED with Jay (2026-07-10):
  each recorded MD step is 25 atomic units of time apart, so
  ts_to_fs = 25 * 0.0241888 fs/au = 0.60472 fs/step (the inline comment
  literally saying "1 au of time in femtoseconds" was imprecise wording
  for this - the constant 6.04721e-16 s already bakes in the factor of
  25 - but the numeric value itself is correct).

Assumptions / non-obvious logic
--------------------------------
- Theta branch correction (theta -> 180 - theta when raw theta > 100
  degrees): CONFIRMED with Jay (2026-07-10) to be an intentional,
  empirical fix. theta/phi are computed upstream (in the trajectory-
  generation step, not in this file) as
  `theta_deg = np.degrees(np.arccos(cos_theta))`, whose principal value
  range is inherently [0, 180] but which can jump discontinuously
  between physically-adjacent frames - e.g. from a sign/orientation
  ambiguity in whatever reference axis cos_theta is measured against, or
  numerical sensitivity of arccos near cos_theta = +/-1 (where its
  derivative diverges). The `> 100` fold empirically restores a
  continuous trajectory for this dataset; it is not a claim about a
  symmetry of the (theta, phi) energy surface itself. Left unchanged.
- Bilinear interpolation via scipy's RegularGridInterpolator (rather than
  the nearest-neighbor KDTree lookup used in plot_para_meta_trajectory_overlay.py) is used to
  evaluate the Para/Meta surface at each MD trajectory point. This is a
  smoother, more accurate mapping than nearest-neighbor, but note
  `fill_value=None` makes the interpolator *extrapolate* (not error or
  NaN) for any (theta, phi) trajectory point that falls outside the
  surface grid's sampled range - such points would silently receive
  linearly-extrapolated (and potentially unphysical) energies with no
  warning. Worth checking trajectory (theta, phi) ranges against the
  grid bounds before trusting results near the domain edges.
- Shading/threshold: +/-5 kcal/mol is hardcoded in
  plot_deltaE_timeseries as the stabilization/destabilization cutoff,
  matching the project's stated classification scheme. (Originally
  hardcoded as +/-6 kcal/mol; corrected to +/-5 per Jay's confirmation
  2026-07-10.)
- An early, unused, commented-out version of main() (using the KDTree
  nearest-neighbor approach) is left below the active implementation;
  kept as-is/untouched since removing it would be a structural change
  beyond cosmetic annotation, but noted here for clarity.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import sys

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_B_wb97x_d_4000_ts.xyz" #"md_trajectory.txt"  # fixed 2026-07-10: was missing direction letter
ENERGY_FILE = "isomer_Nel_49_Nph_10_total_energies.dat"  # fixed 2026-07-10: was a nonexistent placeholder path
AU_TO_KCAL = 627.509  # Hartree -> kcal/mol conversion factor

from scipy.interpolate import RegularGridInterpolator

def get_interpolated_energies(df_grid, df_md):
    """
    Interpolates QED-CCSD energies onto MD trajectory points.

    Builds a bilinear RegularGridInterpolator over the (theta, phi) grid
    for Para_E and Meta_E, evaluates it at every (theta, phi) pair in the
    MD trajectory, and stores the result plus the Para-Meta difference
    (in kcal/mol) on df_md.
    """
    # 1. Reshape the coarse grid into regular 2D arrays
    # Assuming your grid is theta (num_t) x phi (num_p)
    theta_vals = np.sort(df_grid['theta'].unique())
    phi_vals = np.sort(df_grid['phi'].unique())

    # Create the interpolator function for Para and Meta
    # 'linear' is usually best for energy surfaces to avoid artificial oscillations
    # NOTE: fill_value=None -> out-of-bounds (theta, phi) points are
    # linearly extrapolated rather than raising an error or returning NaN.
    interp_para = RegularGridInterpolator((theta_vals, phi_vals),
                                           df_grid.pivot(index='theta', columns='phi', values='Para_E').values,
                                           bounds_error=False, fill_value=None)

    interp_meta = RegularGridInterpolator((theta_vals, phi_vals),
                                          df_grid.pivot(index='theta', columns='phi', values='Meta_E').values,
                                          bounds_error=False, fill_value=None)

    # 2. Extract MD coordinates as a list of (theta, phi) pairs
    pts = df_md[['theta', 'phi']].values

    # 3. Evaluate the interpolator at these points
    df_md['para_e_interp'] = interp_para(pts)
    df_md['meta_e_interp'] = interp_meta(pts)

    # 4. Calculate smoothed difference (Hartree -> kcal/mol)
    df_md['diff_om_grid'] = (df_md['para_e_interp'] - df_md['meta_e_interp']) * AU_TO_KCAL

    return df_md


def parse_md_data(filename):
    """Parses MD block format to extract Step, Energy, Phi, and Theta.

    Reads each "Step ..." header line (e.g.
    "Step 0 E=-436... phi=36.879 theta=73.569") and applies the theta
    branch correction (theta -> 180 - theta when raw theta > 100 deg;
    see module docstring) before returning a per-frame DataFrame.
    """
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

def plot_deltaE_timeseries(df_md, file_name="deltaE_vs_time_direction_B.png"):
    """Publication-ready plot with threshold shading.

    Plots DeltaE_para-meta(t) in kcal/mol together with the trajectory's
    relative electronic energy (offset by its running minimum), shading
    the region where Para is stabilized (<= -5 kcal/mol) or destabilized
    (>= +5 kcal/mol) relative to Meta. The first 10 recorded frames are
    dropped from the plotted window (equilibration/startup transient).
    """
    # 1. Styling
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "font.size": 16,
        "mathtext.fontset": "stix",
    })

    fig, ax = plt.subplots(figsize=(7.5, 5))
    # convert between time step number and femtoseconds
    # CONFIRMED (2026-07-10): each recorded MD step is 25 atomic units of
    # time apart, so ts_to_fs = 25 * 0.0241888 fs/au = 0.60472 fs/step.
    ts_to_fs = 6.04721e-16 * 1e15 # 25 au/step -> fs/step
    # Data extraction (skip first 10 frames - startup/equilibration window)
    t_values = df_md['step'][10:] * ts_to_fs # Example: first 100 steps
    traj1 = df_md['diff_om_grid'][10:] # Para-Meta
    traj2 = (df_md['e_md'] - df_md['e_md'].min())[10:] * AU_TO_KCAL # Rel Energy

    # maroon hex color: #800000
    color1, color2, color3 = "#005F73", "#BB3E03", "#800000"

    # 2. Shading for stabilization/destabilization
    # +/-5 kcal/mol classification threshold (corrected 2026-07-10 from an
    # originally-hardcoded +/-6 to match the project's stated definition).
    ax.fill_between(t_values, traj1, -5, where=(traj1 <= -5), color=color1, alpha=0.15, interpolate=True)
    ax.fill_between(t_values, traj1, 5, where=(traj1 >= 5), color=color3, alpha=0.15, interpolate=True)

    # 3. Plots
    #ax.plot(t_values, traj1, color=color1, linewidth=2.0, alpha=0.9, label="Para vs Meta ($\Delta E$)")
    ax.plot(t_values, traj2, color=color2, linewidth=2.0, alpha=0.9, label="Nitrobenzene Relative Energy (kcal / mol)")
    ax.plot(t_values, traj1, color=color1, linewidth=2.0, alpha=0.9, label="$\Delta E_{para-meta}$ (kcal / mol)")
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

# --- Earlier draft of main(), kept for reference (uses KDTree nearest-neighbor
# mapping instead of interpolation). Left commented out/unused, unchanged
# from the original script. ---
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
    """Load MD trajectory + energy grid, interpolate, and plot the time series."""
    # 1. Load and Clean
    df_md = parse_md_data(MD_FILE)
    df_grid = pd.read_csv(ENERGY_FILE, sep='\s+', skiprows=[1])

    # Force numeric for grid
    for col in ['theta', 'phi', 'Para_E', 'Meta_E']:
        df_grid[col] = pd.to_numeric(df_grid[col], errors='coerce')

    # 2. Interpolate Energies (Replaces the KDTree block)
    df_md = get_interpolated_energies(df_grid, df_md)

    # 3. Generate Plot
    plot_deltaE_timeseries(df_md)

if __name__ == "__main__":
    main()
