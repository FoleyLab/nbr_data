# ΔE Timeseries & Basin Dwell-Time Analysis Along AIMD Trajectories

**Data for:** Figure 6 — Molecular Dynamics Relaxed Timeseries

**System:** Nitrobenzene (14 atoms) cavity-coupled AIMD trajectories (ωB97X-D, 4000 steps) along field **direction A** (θ = 70°, φ = 31°) and **direction D** (θ = 65°, φ = 78°). The relaxed QED-CCSD energy differences ΔE(Ortho − Meta) and ΔE(Para − Meta) from Figure 2 are interpolated onto each trajectory's (θ, φ) path, and every frame is classified into a stabilization basin to quantify how long the molecule dwells in each cavity-favored regime.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Directory Overview](#directory-overview)
3. [Analysis Details](#analysis-details)
4. [Reproducibility & Plotting](#reproducibility--plotting)
5. [Swapping in New Intermediate-Energy Data](#swapping-in-new-intermediate-energy-data)
6. [Citation](#citation)

---

## Quick Start

| File | Description |
|------|-------------|
| `deltaE_vs_time_direction_A.png` | ΔE vs. time with basin shading, direction A |
| `deltaE_vs_time_direction_D.png` | ΔE vs. time with basin shading, direction D |

To regenerate the figure **and** the console dwell-time report:

```bash
python plot_timeseries_with_dwell_times.py
```

The script is currently configured for **direction D**; edit `MD_FILE` and the output filename for direction A (see [Reproducibility & Plotting](#reproducibility--plotting)).

---

## Directory Overview

```
Figure_06_Molecular_Dynamics_Relaxed_Timeseries/
├── README.md                                       # This file
├── nitrobenzene_direction_A_wb97x_d_4000_ts.xyz    # AIMD trajectory, direction A (4000 steps)
├── nitrobenzene_direction_D_wb97x_d_4000_ts.xyz    # AIMD trajectory, direction D (4000 steps)
├── Relaxed_QED_CCSD_Combined_Results.txt           # Energy-surface grid (copy from Figure_02)
├── plot_timeseries_with_dwell_times.py             # Script: timeseries plot + dwell-time report
└── deltaE_vs_time_direction_{A,D}.png              # Output figures
```

The trajectory and energy-grid file formats are identical to those documented in `Figure_05_Molecular_Dynamics_on_Relaxed_Surface/README.md` (multi-frame XYZ with `Step <n> E=… phi=… theta=…` comment lines; whitespace-delimited θ-major grid with columns `theta phi Ex Ey Ez Para_E Ortho_E Meta_E`, Hartree).

---

## Analysis Details

### Time axis

Step index → femtoseconds via `TS_TO_FS`: 25 au/step ≈ **0.604721 fs/step** (confirmed in the script's module docstring).

### Surface interpolation

Trajectory (θ, φ) points are mapped onto the energy grid by bilinear interpolation (`RegularGridInterpolator`, `fill_value=None`), which extrapolates rather than erroring for points outside the sampled grid range.

### Theta branch correction

Raw θ values above 100° are mapped as θ → 180° − θ. This is an intentional, empirical fix for arccos branch jumps in the upstream angle calculation (`np.arccos` principal range vs. reference-axis sign ambiguity), confirmed 2026-07-10; it is not a statement about surface symmetry.

### Basin classification & dwell times

Each frame is classified from the two interpolated differences (in kcal/mol), using threshold `THRESH`:

| Basin | Criterion |
|-------|-----------|
| **Ortho-favored** | ΔE(ortho−meta) ≤ −THRESH **and** ΔE(ortho−meta) < ΔE(para−meta) |
| **Para-favored** | ΔE(para−meta) ≤ −THRESH **and** ΔE(para−meta) < ΔE(ortho−meta) |
| **Meta-favored** | ΔE(ortho−meta) ≥ +THRESH **and** ΔE(para−meta) ≥ +THRESH |
| **Neutral** | none of the above |

Consecutive frames in the same state are grouped into contiguous runs (`itertools.groupby`); each run's dwell time is (frame count) × (fs/step). The chronological printout lists Ortho and Para dwell periods (Neutral runs are computed but not printed). The figure shades winner-take-all regions (Ortho-wins / Para-wins / both-destabilized) with the same criteria.

> **Threshold caveat:** the classification threshold is the module constant `THRESH` (kcal/mol). The module docstring describes the project's ±5 kcal/mol classification scheme (corrected from an earlier ±6 hardcoding, per confirmation on 2026-07-10), while the constant in the current code is set to `1`. Check `THRESH` before regenerating dwell-time numbers, and keep figures/tables consistent with the value stated in the manuscript.

---

## Reproducibility & Plotting

```bash
python plot_timeseries_with_dwell_times.py
```

Configuration block (top of script):

| Setting | Current value | Notes |
|---------|---------------|-------|
| `MD_FILE` | `nitrobenzene_direction_D_wb97x_d_4000_ts.xyz` | Swap to the direction-A file for that figure |
| `ENERGY_FILE` | `Relaxed_QED_CCSD_Combined_Results.txt` | See next section (placeholder path fixed 2026-07-10) |
| `AU_TO_KCAL` | 627.509 | Hartree → kcal/mol |
| `TS_TO_FS` | 0.604721 | 25 au/step → fs/step |
| `THRESH` | 1 | Basin threshold, kcal/mol — see caveat above |
| output | `deltaE_vs_time_direction_D.png` (600 dpi) | Update suffix when switching directions |

Console output: chronological list of dwell times (fs) in the Ortho- and Para-favored basins.

---

## Swapping in New Intermediate-Energy Data

To run the analysis on an energy surface from a different level of theory (e.g. the pQED `isomer_Nel_49_Nph_10_total_energies.dat` grid used elsewhere in this project), point `ENERGY_FILE` at the new file. It must be whitespace-delimited with columns `theta phi ... Para_E Ortho_E Meta_E` (Hartree), two header lines, and rows in θ-major order (see the Figure_05 README for details).

---

## Citation

> *Cite the manuscript (DOI to be added upon publication). This dataset is hosted on Zenodo as supporting data for the associated publication.*

**Keywords:** ab initio molecular dynamics · cavity QED · dwell time · energy-difference timeseries · nitrobenzene · Wheland intermediate
