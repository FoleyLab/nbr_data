"""
Trajectory-over-PES overlay plots: BOTH Ortho-Meta and Para-Meta surfaces.

Purpose
-------
Despite the "_para_meta" filename, this script generates TWO overlay
figures per trajectory: one on the Ortho_E - Meta_E surface and one on
the Para_E - Meta_E surface (see flagged issues re: misleading name).
Each figure shows the static energy-difference landscape as a heatmap
with the AIMD trajectory's (theta, phi) path drawn over it, plus
start/end markers.

Inputs
------
- MD_FILE: AIMD trajectory in the project's ".xyz + Step header" format
  (currently hardcoded to direction C). Each frame's header line looks
  like "Step <n>  E=<hartree>  phi=<deg>  theta=<deg>".
- ENERGY_FILE: whitespace-delimited surface grid with columns
      theta  phi  Ex  Ey  Ez  Para_E  Ortho_E  Meta_E
  (Hartree). RESOLVED: was set to a nonexistent "CCSD_Combined_Results.txt";
  confirmed with Jay 2026-07-10 and set to
  "isomer_Nel_49_Nph_10_total_energies.dat" (90x90 grid, Nel=49
  electronic states / Nph=10 photon Fock states). See README's
  "Swapping in new intermediate-energy data" section to point this at a
  different level of theory later.

Outputs
-------
- traj_overlay_ortho_meta_direction_C.png
- traj_overlay_para_meta_direction_C.png

Units
-----
- Angles in degrees; energies read in Hartree and converted to kcal/mol
  (AU_TO_KCAL) only for the plotted quantities.

Assumptions / non-obvious logic
--------------------------------
- Theta branch correction (theta -> 180 - theta when raw theta > 100
  degrees): CONFIRMED with Jay (2026-07-10) as an intentional, empirical
  fix for arccos branch jumps in the upstream angle calculation (theta
  is computed there as np.degrees(np.arccos(cos_theta)), whose principal
  range is [0, 180] but which can jump discontinuously between frames
  due to reference-axis sign ambiguity or numerical sensitivity near
  cos_theta = +/-1). Not a statement about surface symmetry. Left
  unchanged here; see plot_para_meta_trajectory_overlay.py's module docstring for the full
  explanation.
- Grid reshape for the background surfaces assumes df_grid's rows are
  already ordered as a perfect theta-major grid (num_t x num_p). This is
  more fragile than the pivot-table approach used in plot_para_meta_timeseries.py /
  plot_ortho_para_timeseries.py / plot_timeseries_with_dwell_times.py, which sorts explicitly and
  would tolerate a differently-ordered input file; flagged as an
  inconsistency across scripts.
- Color scale is symmetric about zero (vlimit = max abs value) so that
  stabilizing (negative, cool) and destabilizing (positive, warm) shifts
  are visually comparable across both surfaces.
- The MD trajectory line itself is drawn directly from the parsed
  (phi, theta) pairs - it does not need interpolated energies, only
  angular position, since it is overlaid on top of the surface rather
  than plotted as its own energy trace.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_C_wb97x_d_4000_ts.xyz"
ENERGY_FILE = "isomer_Nel_49_Nph_10_total_energies.dat"  # fixed 2026-07-10: was a nonexistent placeholder path
AU_TO_KCAL = 627.509  # Hartree -> kcal/mol conversion factor
DPI = 600 # High resolution for publication

def parse_md_data(filename):
    """Parses MD block format and applies theta branch correction.

    Extracts (step, phi, theta) from each "Step ..." header line. Theta
    is folded via theta -> 180 - theta whenever the raw parsed value
    exceeds 100 degrees (confirmed, intentional arccos branch-jump
    correction; see module docstring).
    """
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
    """Builds and saves the Ortho-Meta and Para-Meta trajectory-overlay figures."""
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
    # NOTE: assumes df_grid rows form a perfect theta-major grid of shape
    # (num_t, num_p); see module docstring caveat.
    unique_t = np.unique(df_grid['theta'])
    unique_p = np.unique(df_grid['phi'])
    num_t, num_p = len(unique_t), len(unique_p)

    T = df_grid['theta'].values.reshape(num_t, num_p)
    P = df_grid['phi'].values.reshape(num_t, num_p)

    # Calculate Difference Surfaces (Hartree -> kcal/mol via AU_TO_KCAL)
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
