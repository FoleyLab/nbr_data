"""
Time series plot + basin dwell-time analysis, along a single AIMD
trajectory.

Purpose
-------
The most complete of the timeseries scripts: interpolates Ortho/Para/
Meta energies onto an AIMD trajectory, classifies every frame into an
Ortho-preferred / Para-preferred / Meta-preferred / Neutral basin per the
project's stabilization criteria, prints a chronological breakdown of how
long each visit to a basin lasted ("dwell time"), and produces the same
style of DeltaE-vs-time figure as plot_ortho_para_timeseries.py with comparative
(winner-take-all) shading.

Inputs
------
- MD_FILE: AIMD trajectory, ".xyz + Step header" format (hardcoded to
  direction D here).
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
- Console printout: chronological list of dwell times (fs) in each basin
  (Ortho/Para only; Neutral periods are computed but not printed).
- deltaE_vs_time_direction_D.png

Units
-----
- Angles in degrees; energies converted Hartree -> kcal/mol (AU_TO_KCAL).
- Time axis: step index -> femtoseconds via TS_TO_FS / ts_to_fs.
  CONFIRMED with Jay (2026-07-10): each recorded MD step is 25 atomic
  units of time apart, so TS_TO_FS = 25 * 0.0241888 fs/au = 0.60472
  fs/step.

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
- Bilinear interpolation (RegularGridInterpolator, fill_value=None)
  extrapolates rather than erroring for MD (theta, phi) points outside
  the surface grid's sampled range - see plot_para_meta_timeseries.py's docstring.
- Basin classification / dwell-time bookkeeping (calculate_dwell_times):
    * Ortho-favored:  DeltaE_ortho-meta <= -5  AND  DeltaE_ortho-meta < DeltaE_para-meta
    * Para-favored:   DeltaE_para-meta  <= -5  AND  DeltaE_para-meta  < DeltaE_ortho-meta
    * Meta-favored:   DeltaE_ortho-meta >= +5  AND  DeltaE_para-meta  >= +5 (both destabilized)
    * Neutral:        none of the above
  Matches the project's stated +/-5 kcal/mol classification threshold.
  (Originally hardcoded as +/-6 kcal/mol in four places; corrected to
  +/-5 per Jay's confirmation 2026-07-10 - this changes which frames are
  classified into which basin and will shift the printed dwell times
  relative to any previously-generated output using the old threshold.)
  Frame-to-frame states are grouped into contiguous runs with
  itertools.groupby (so consecutive frames of the same state count as
  one "visit"/dwell period), and each run's duration is
  (frame count) * (fs per step). Only Ortho/Para dwell periods are
  printed (Neutral runs are computed but skipped in the printout, per
  the original script).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from itertools import groupby
import sys

# --- Configuration ---
MD_FILE = "nitrobenzene_direction_D_wb97x_d_4000_ts.xyz"
ENERGY_FILE = "isomer_Nel_49_Nph_10_total_energies.dat"  # fixed 2026-07-10: was a nonexistent placeholder path
AU_TO_KCAL = 627.509  # Hartree -> kcal/mol conversion factor
TS_TO_FS = 6.04721e-16 * 1e15 # 25 au/step -> fs/step (confirmed, see module docstring)
THRESH = 1 # threshold in kcal / mol for stabilization


def calculate_dwell_times(df_md, ts_to_fs):
    """Identifies and prints the duration spent in each stabilizing basin.

    Classifies every trajectory frame into "Ortho", "Para", "Meta", or
    "Neutral" using the +/-5 kcal/mol comparative thresholds (matches the
    project's stated classification scheme; corrected 2026-07-10 from an
    originally-hardcoded +/-6), then groups consecutive same-state frames
    into discrete dwell periods and prints their duration in
    femtoseconds, in chronological order.
    """

    # 1. Assign a state to every timestep based on your logic
    # Using 'N' for Neutral (between -5 and 5)
    states = []
    traj_om = df_md['diff_om'].values
    traj_pm = df_md['diff_pm'].values
    
    for om, pm in zip(traj_om, traj_pm):
        if om <= -THRESH and om < pm:
            states.append("Ortho")
        elif pm <= -THRESH and pm < om:
            states.append("Para")
        elif om >= THRESH and pm >= THRESH:
            states.append("Meta")
        else:
            states.append("Neutral")

    # 2. Group contiguous states and calculate duration
    # groupby returns (key, group_iterator); consecutive frames sharing
    # the same state are treated as a single continuous "dwell" period.
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
    """Interpolates QED-CCSD energies (Ortho, Para, Meta) onto MD trajectory.

    Builds one bilinear RegularGridInterpolator per intermediate over the
    (theta, phi) surface grid and evaluates each at every MD trajectory
    point, storing the raw interpolated energies and the two energy
    differences of interest (Ortho-Meta, Para-Meta) in kcal/mol.
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
    100 degrees (see module docstring caveat).
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

def plot_deltaE_timeseries(df_md, file_name="deltaE_vs_time_direction_D.png"):
    """Plot with comparative shading to show the thermodynamic winner.

    Plots DeltaE_ortho-meta(t), DeltaE_para-meta(t), and relative
    electronic energy, shading regions by which basin is favored
    (Ortho-preferred, Para-preferred, or both intermediates
    destabilized/"Meta preferred") using the same +/-5 kcal/mol
    comparative logic as calculate_dwell_times. The first 10 recorded
    frames are dropped (equilibration/startup transient).
    """
    # 1. Styling
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "font.size": 16,
        "mathtext.fontset": "stix",
    })

    fig, ax = plt.subplots(figsize=(10, 6))

    # Time conversion: 25 au/step -> fs/step (confirmed, see module docstring)
    ts_to_fs = 6.04721e-16 * 1e15

    # Data Windowing (skip first 10 frames - startup/equilibration window)
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
    # +/-5 kcal/mol classification threshold, matching calculate_dwell_times
    # (corrected 2026-07-10 from an originally-hardcoded +/-6 to match the
    # project's stated definition).

    # ORTHO WINNER: Ortho < -5 AND Ortho < Para
    ortho_wins = (traj_om <= -THRESH) & (traj_om < traj_pm)
    ax.fill_between(t_fs, traj_om, -THRESH, where=ortho_wins,
                    color=c_ortho, alpha=0.2, interpolate=True, label='Ortho Preferred')

    # PARA WINNER: Para < -5 AND Para < Ortho
    para_wins = (traj_pm <= -THRESH) & (traj_pm < traj_om)
    ax.fill_between(t_fs, traj_pm, -THRESH, where=para_wins,
                    color=c_para, alpha=0.2, interpolate=True, label='Para Preferred')

    # BOTH DESTABILIZED: Both > 5
    both_bad = (traj_om >= THRESH) & (traj_pm >= THRESH)
    ax.fill_between(t_fs, THRESH, np.maximum(traj_om, traj_pm), where=both_bad,
                    color=c_destab, alpha=0.15, interpolate=True, label='Meta Preferred')

    # 3. Trajectory Lines
    ax.plot(t_fs, traj_rel_e, color=c_nitro, lw=1.5, alpha=0.7, label="Nitrobenzene Rel. $E$")
    ax.plot(t_fs, traj_om, color=c_ortho, lw=2.2, label=r"$\Delta E_{ortho-meta}$")
    ax.plot(t_fs, traj_pm, color=c_para, lw=2.2, linestyle='--', label=r"$\Delta E_{para-meta}$")

    # 4. Ref lines & Formatting
    ax.axhline(0.0, color='black', linestyle='-', linewidth=1.2, alpha=0.7, zorder=0)
    ax.axhline(5.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.axhline(-5.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)

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
    """Load MD trajectory + energy grid, interpolate, report dwell times, and plot."""
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
