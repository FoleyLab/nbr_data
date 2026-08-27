# Cavity vs. Cavity-Free Potential-Energy Fluctuations Along an AIMD Trajectory

**Data for:** Figure S8 — Nitrobenzene Energy Fluctuations

**System:** Nitrobenzene (14 atoms) sampled along the direction-A AIMD trajectory (field orientation θ = 70°, φ = 31°; see `Figure_05_Molecular_Dynamics_on_Relaxed_Surface/`). Cavity-coupled and cavity-free single-point energies are evaluated at the **same** sampled nuclear geometries Rₜ, so their difference isolates the instantaneous (and time-averaged) energetic shift induced by the cavity.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Directory Overview](#directory-overview)
3. [Data File Format Reference](#data-file-format-reference)
4. [How the Data Were Generated](#how-the-data-were-generated)
5. [Citation](#citation)

---

## Quick Start

| File | Description |
|------|-------------|
| `cavity_vs_no_cavity_potential_energy.png` / `.pdf` | Final figure: cavity vs. cavity-free potential energy along the trajectory |

To regenerate the figure from the CSV:

```bash
python plot_cavity_vs_no_cavity_potential.py
```

---

## Directory Overview

```
Figure_S8_Nitrobenzene_Energy_Fluctuations/
├── README.md                                                # This file
├── nitrobenzene_direction_A_wb97x_d_cavity_free_energies.csv# Primary data (1001 sampled frames)
├── cqed_dft_energy_gradient.py                              # CQED-DFT energy/gradient template (cavity on)
├── dft_energy_template.py                                   # Driver: cavity-free recomputation → CSV
├── plot_cavity_vs_no_cavity_potential.py                    # Plotting script
└── cavity_vs_no_cavity_potential_energy.{png,pdf}           # Output figure
```

---

## Data File Format Reference

### `nitrobenzene_direction_A_wb97x_d_cavity_free_energies.csv`

One row per sampled trajectory frame (1001 rows + header). Frames are sampled from the 4000-step direction-A trajectory; the `frame` index maps to MD step = `frame` − 1 (row 0 carries the initial direction-A reference frame, `Step -1 ... phi=31.0 theta=70.0`).

| Column | Description |
|--------|-------------|
| `frame` | Sample index (MD step = frame − 1) |
| `cavity_energy_hartree` | Cavity-coupled energy at Rₜ (Hartree) |
| `no_cavity_energy_hartree` | Cavity-free energy at the same Rₜ (Hartree) |
| `delta_cavity_minus_no_cavity_hartree` | Difference (Hartree) |
| `delta_cavity_minus_no_cavity_kcal_mol` | Difference (kcal/mol; 1 Hartree = 627.509 kcal/mol) |
| `xyz_comment` | Source frame's comment line (`Step <n> E=… phi=… theta=…`) |

---

## How the Data Were Generated

| Script | Role |
|--------|------|
| `cqed_dft_energy_gradient.py` | Template for the cavity side: CQED-DFT energy and nuclear gradient via `cqed_scf`'s `CQEDCalculator` / `CQEDConfig` (Psi4), with an embedded nitrobenzene geometry |
| `dft_energy_template.py` | Driver for the cavity-free side: parses each frame's geometry **and** the cavity energy already stored in the trajectory XYZ comment lines, recomputes the same geometry with zero cavity field using Psi4 (parallelized with `ProcessPoolExecutor`), and writes both energies plus their difference to the CSV |
| `plot_cavity_vs_no_cavity_potential.py` | Reads the CSV and produces the SI figure (PNG + PDF) |

The source trajectory (`nitrobenzene_direction_A_wb97x_d_4000_ts.xyz`, ωB97X-D) is stored in `Figure_05_Molecular_Dynamics_on_Relaxed_Surface/`.

---

## Citation

> *Cite the manuscript (DOI to be added upon publication). This dataset is hosted on Zenodo as supporting data for the associated publication.*

**Keywords:** cavity QED · QED-DFT · potential-energy fluctuations · cavity-free reference · nitrobenzene · molecular dynamics
