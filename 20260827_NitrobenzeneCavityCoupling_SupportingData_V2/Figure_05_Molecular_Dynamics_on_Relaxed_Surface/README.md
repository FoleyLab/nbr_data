# AIMD Trajectories Overlaid on the Relaxed QED-CCSD Energy Surfaces

**Data for:** Figure 5 — Molecular Dynamics on the Relaxed Surface

**System:** Nitrobenzene (14 atoms) propagated by cavity-coupled AIMD (ωB97X-D) with the cavity field oriented along **direction A** (θ = 70°, φ = 31°) or **direction D** (θ = 65°, φ = 78°). Each trajectory's instantaneous cavity-orientation angles (θ, φ) are overlaid on the relaxed QED-CCSD energy-difference landscapes (Ortho − Meta and Para − Meta) from Figure 2.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Directory Overview](#directory-overview)
3. [Data File Format Reference](#data-file-format-reference)
4. [Reproducibility & Plotting](#reproducibility--plotting)
5. [Swapping in New Intermediate-Energy Data](#swapping-in-new-intermediate-energy-data)
6. [Citation](#citation)

---

## Quick Start

| File | Description |
|------|-------------|
| `traj_overlay_ortho_meta_direction_A.png` | Direction-A trajectory over the ΔE (Ortho − Meta) surface |
| `traj_overlay_para_meta_direction_A.png` | Direction-A trajectory over the ΔE (Para − Meta) surface |
| `traj_overlay_ortho_meta_direction_D.png` | Direction-D trajectory over the ΔE (Ortho − Meta) surface |
| `traj_overlay_para_meta_direction_D.png` | Direction-D trajectory over the ΔE (Para − Meta) surface |

To regenerate the plots:

```bash
python plot_ortho_para_trajectory_overlays.py
```

The script is currently configured for **direction D**. To produce the direction-A overlays, edit `MD_FILE` and the two output filename stems in the configuration block (see [Reproducibility & Plotting](#reproducibility--plotting)).

---

## Directory Overview

```
Figure_05_Molecular_Dynamics_on_Relaxed_Surface/
├── README.md                                       # This file
├── nitrobenzene_direction_A_wb97x_d_4000_ts.xyz    # AIMD trajectory, direction A (4000 steps)
├── nitrobenzene_direction_D_wb97x_d_4000_ts.xyz    # AIMD trajectory, direction D (4000 steps)
├── Relaxed_QED_CCSD_Combined_Results.txt           # Energy-surface grid (copy from Figure_02)
├── plot_ortho_para_trajectory_overlays.py          # Script: overlay plots for both surfaces
├── traj_overlay_ortho_meta_direction_{A,D}.png     # Output: Ortho − Meta overlays
└── traj_overlay_para_meta_direction_{A,D}.png      # Output: Para − Meta overlays
```

Despite the script's filename (`..._para_meta`), it generates **both** the Ortho − Meta and Para − Meta overlays in a single run.

---

## Data File Format Reference

### Trajectory files (`nitrobenzene_direction_*_wb97x_d_4000_ts.xyz`)

Standard multi-frame XYZ: 14 atoms per frame, 4000 frames. Each frame's comment line carries the step index, the cavity-coupled energy, and the instantaneous field-orientation angles of the NO₂ group:

```
14
Step <n>  E=<hartree>  phi=<deg>  theta=<deg>
C  x  y  z
...
```

### Energy-surface file (`Relaxed_QED_CCSD_Combined_Results.txt`)

Copy of the symmetry-expanded relaxed QED-CCSD grid from `Figure_02_Relaxed_Theta_Phi_Scan/`. Whitespace-delimited with **two header lines**; columns:

| Column | Description |
|--------|-------------|
| `theta` | Polar angle (degrees, 0–180) |
| `phi` | Azimuthal angle (degrees, 0–360) |
| `Ex`, `Ey`, `Ez` | Field unit-vector components |
| `Para_E` | CCSD total energy, para (Hartree) |
| `Ortho_E` | CCSD total energy, ortho (Hartree) |
| `Meta_E` | CCSD total energy, meta (Hartree) |

Rows are ordered as a perfect θ-major grid (num_θ × num_φ); the script reshapes columns directly onto a (num_θ, num_φ) mesh. Energies are converted to kcal/mol with 1 Hartree = 627.509 kcal/mol.

---

## Reproducibility & Plotting

```bash
python plot_ortho_para_trajectory_overlays.py
```

Configuration block (top of script):

| Setting | Current value | Notes |
|---------|---------------|-------|
| `MD_FILE` | `nitrobenzene_direction_D_wb97x_d_4000_ts.xyz` | Swap to the direction-A file for the A overlays |
| `ENERGY_FILE` | `Relaxed_QED_CCSD_Combined_Results.txt` | See next section |
| output stems | `traj_overlay_{ortho,para}_meta_direction_D.png` | Update the suffix when switching directions |

The direction-A PNGs in this folder were produced by an earlier run with `MD_FILE` and the output stems set accordingly.

> **Known caveat:** the script's module docstring is stale in places — it refers to "direction C" (a holdover from an earlier configuration) and to a since-replaced placeholder energy file. The configuration block above reflects the actual current state (confirmed against the code, 2026-07-10).

---

## Swapping in New Intermediate-Energy Data

To overlay the trajectories on an energy surface from a different level of theory (e.g. the pQED `isomer_Nel_49_Nph_10_total_energies.dat` grid used elsewhere in this project):

1. Point `ENERGY_FILE` at the new file.
2. Ensure it provides the columns `theta phi ... Para_E Ortho_E Meta_E` (extra columns such as `Ex Ey Ez` are fine), whitespace-delimited, with **two header lines** (the reader uses `skiprows=[1]` after the column header) and rows in **θ-major order** — the script reshapes the energy columns directly onto a (num_θ, num_φ) grid without re-sorting.
3. Energies are assumed to be in Hartree and are converted to kcal/mol internally.

---

## Citation

> *Cite the manuscript (DOI to be added upon publication). This dataset is hosted on Zenodo as supporting data for the associated publication.*

**Keywords:** ab initio molecular dynamics · cavity QED · energy-difference landscape · trajectory overlay · nitrobenzene · Wheland intermediate
