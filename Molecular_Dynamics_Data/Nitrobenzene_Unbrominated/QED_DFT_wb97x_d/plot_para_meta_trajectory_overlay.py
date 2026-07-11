"""
Trajectory-over-PES overlay plot: Para vs. Meta energy-difference surface.

Purpose
-------
Overlays a single AIMD trajectory's (theta, phi) path onto the static
Para_E - Meta_E energy-difference surface (the "landscape" the cavity
orientation angles are sampled over), producing a 2D heatmap with the
trajectory drawn on top (start/end markers included).

Inputs
------
- MD_FILE: an AIMD trajectory file in the project's ".xyz + Step header"
  format. Each frame is preceded by a line such as
      "Step <n>  E=<hartree>  phi=<deg>  theta=<deg>"
  from which (step, energy, theta, phi) are parsed.
- ENERGY_FILE: whitespace-delimited surface data file with a header row,
  a divider row (hence skiprows=[1]), and columns
      theta  phi  Ex  Ey  Ez  Para_E  Ortho_E  Meta_E
  (energies in Hartree). This is the (theta, phi) grid the Pauli-Fierz /
  EOM-CCSD (or replacement level of theory) ground-state energies were
  computed on for each Wheland intermediate (see project README).

  RESOLVED: was pointing at a nonexistent "CCSD_Combined_Results.txt";
  confirmed with Jay 2026-07-10 and set to
  "isomer_Nel_49_Nph_10_total_energies.dat" (90x90 grid, Nel=49
  electronic states / Nph=10 photon Fock states). See README's
  "Swapping in new intermediate-energy data" section to point this at a
  different level of theory later.

Output
------
- traj_overlay_direction_D.png: PNG figure of the trajectory over the
  Para-Meta surface for whichever direction MD_FILE points at (currently
  hardcoded to direction D; see flagged issues re: making this a
  parameter instead of editing the script per run).

Units
-----
- Angles (theta, phi) in degrees throughout.
- Energies read in Hartree; converted to kcal/mol via AU_TO_KCAL for
  plotting only (the underlying .dat/.txt files remain in Hartree).

Assumptions / non-obvious logic
--------------------------------
- Bilinear interpolation: MD (theta, phi) points are mapped onto the
  Para_E/Meta_E surface via scipy's RegularGridInterpolator. CONSOLIDATED
  2026-07-10 (per Jay's request) to match the interpolation approach used
  in plot_para_meta_timeseries.py / plot_ortho_para_timeseries.py /
  plot_timeseries_with_dwell_times.py, replacing the previous KDTree
  nearest-neighbor lookup so all five scripts now use one consistent
  method. As in those scripts, `fill_value=None` means out-of-grid-bounds
  (theta, phi) points are linearly extrapolated rather than erroring.
  (This computed column is not currently plotted by this script - see
  below - so this change does not alter the saved figure, only the
  method used to compute the otherwise-unused interpolated values.)
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
- The grid reshape (T, P, diff_om below) assumes the rows of ENERGY_FILE
  are already ordered as a perfect theta-major, phi-minor grid matching
  `num_t x num_p`. If that ordering assumption is violated, the reshape
  will silently produce a corrupted-looking surface with no error.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
from scipy.interpolate import RegularGridInterpolator

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_D_wb97x_d_4000_ts.xyz"  # Your MD output file
ENERGY_FILE = "isomer_Nel_49_Nph_10_total_energies.dat"  # fixed 2026-07-10: was a nonexistent placeholder path
AU_TO_KCAL = 627.509  # Hartree -> kcal/mol conversion factor, used for plotting/display only

def parse_md_data(filename):
    """Parses your specific MD block format.

    Reads every line containing "Step" (e.g.
    "Step 0  E=-436.4856647627  phi=36.869  theta=73.569") and extracts
    the step index, total electronic energy (Hartree), and the two
    cavity-orientation angles theta/phi (degrees).

    Applies a theta branch correction: if the raw theta exceeds 100
    degrees, it is folded via theta -> 180 - theta. See module docstring
    "Assumptions" section - this is a confirmed, intentional empirical
    fix for arccos branch jumps in the upstream angle calculation, not
    altered here.
    """
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

# 2. Map MD coordinates onto the Grid via bilinear interpolation
# (consolidated 2026-07-10 to match the RegularGridInterpolator approach
# used in the other four scripts, replacing the previous KDTree
# nearest-neighbor lookup - see module docstring "Assumptions" section.)
theta_vals = np.sort(df_grid['theta'].unique())
phi_vals = np.sort(df_grid['phi'].unique())

def _make_interp(col_name):
    """Builds a bilinear interpolator over the (theta, phi) grid for one energy column."""
    grid_data = df_grid.pivot(index='theta', columns='phi', values=col_name).values
    # fill_value=None -> linearly extrapolates for out-of-bounds (theta, phi)
    # points rather than erroring or returning NaN (see module docstring).
    return RegularGridInterpolator((theta_vals, phi_vals), grid_data,
                                   bounds_error=False, fill_value=None)

interp_para = _make_interp('Para_E')
interp_meta = _make_interp('Meta_E')
pts = df_md[['theta', 'phi']].values

# Get the corresponding Para/Meta energies from the grid
# NOTE: despite the "_om" ("ortho-meta") naming below, this column is actually
# Para_E - Meta_E (see AU_TO_KCAL multiply). It is computed here but not
# subsequently used anywhere in the plot itself (the overlay below only
# needs the trajectory's (phi, theta) positions, not its energy values) -
# left in place unchanged since it doesn't affect the figure produced.
df_md['para_e_grid'] = interp_para(pts)
df_md['meta_e_grid'] = interp_meta(pts)
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
# Reshape the flat (theta, phi, energy) grid rows into 2D arrays suitable
# for pcolormesh. This assumes df_grid rows are already sorted as a
# regular theta-major grid (see module docstring caveat).
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
