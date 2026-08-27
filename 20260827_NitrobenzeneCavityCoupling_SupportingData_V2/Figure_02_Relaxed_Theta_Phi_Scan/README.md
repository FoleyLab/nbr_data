# QED-DFT Geometry Optimizations & QED-CCSD Single-Point Energies on a Spherical Polar Grid

**Data for:** Figure 2 — Theta/Phi Cavity-Orientation Scan of Bromo-Nitrobenzene Wheland Intermediates

**System:** Three tautomeric Wheland intermediates of bromo-nitrobenzene — **ortho**, **meta**, and **para** (relative to the Br substituent) — each (15 atoms, charge +1, doublet-derived singlet) relaxed under constrained QED-DFT at 191 unique cavity-field orientations on a spherical polar (θ, φ) grid, then evaluated at the QED-CCSD level of theory on the relaxed geometries.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Directory Overview](#directory-overview)
3. [The (θ, φ) Grid](#the-grid)
4. [Computational Methods](#computational-methods)
5. [Per-Cell Run Directory Layout](#per-cell-run-directory-layout)
6. [Data File Format Reference](#data-file-format-reference)
7. [Reproducibility & Plotting](#reproducibility--plotting)
8. [Validation & Quality Control](#validation--quality-control)
9. [Citation](#citation)

---

## Quick Start

For collaborators who want the final energy values without digging into the raw calculations:

| File | Description |
|------|-------------|
| `relaxed_qed_ccsd_intermediate_scans.csv` | **Primary data table** — 191 rows (one per unique θ, φ), 8 columns of SCF and CCSD energies for all three isomers |
| `Relaxed_QED_CCSD_Combined_Results.txt` | Symmetry-expanded grid — 440 rows × 8 columns, includes computed field vector components (Eₓ, Eᵧ, E_z) |
| `relaxed_ortho_meta_diff_QED_CCSD_22.png` | Energy difference map: Ortho − Meta (kcal/mol) |
| `relaxed_para_meta_diff_QED_CCSD_22.png` | Energy difference map: Para − Meta (kcal/mol) |

To regenerate the combined results file from the CSV:

```bash
python read_csv_transform_to_full_theta_phi_grid.py
```

To regenerate the plots from the combined results:

```bash
python plot_single_panel_ortho_para_meta.py
```

---

## Directory Overview

```
Figure_02_Revised_Theta_Phi_Scan/
├── README.md                                          # This file
├── relaxed_qed_ccsd_intermediate_scans.csv            # Primary data: 191 rows × 8 energy columns
├── Relaxed_QED_CCSD_Combined_Results.txt              # Symmetry-expanded: 440 rows × 8 columns
├── relaxed_ortho_meta_diff_QED_CCSD_22.png            # Plot: ΔE (Ortho − Meta), kcal/mol
├── relaxed_para_meta_diff_QED_CCSD_22.png             # Plot: ΔE (Para − Meta), kcal/mol
├── read_csv_transform_to_full_theta_phi_grid.py       # Script: expand CSV to full θ-φ grid
├── plot_single_panel_ortho_para_meta.py               # Script: generate Figure 2 plots
│
├── QED_DFT_Theta_Phi_Relax_Campaign/                 # ── QEDEDT relaxation campaign ──
│   ├── relaxed_qed_ccsd_intermediate_scans.csv        # Copy of primary data (QED-CCSD energies on relaxed geoms)
│   ├── unrelaxed_qed_ccsd_intermediate_scans.csv      # Reference: QED-CCSD on pristine (unrelaxed) geoms
│   ├── intermediate_scans.csv                          # The (θ,φ) scan grid (subset)
│   ├── campaign_folder/                                # Campaign scripts and input files
│   │   ├── campaign.py                               # Main optimization driver
│   │   ├── grid_campaign_no_freq.py                  # Grid-scan optimization driver (warm-start chains)
│   │   ├── cavity_common.py                          # Shared QED-DFT physics configuration
│   │   ├── ccsd_input.py                             # Builds QED-CCSD input JSON from relaxed geometry
│   │   ├── campaign_intro.md                         # Documentation for grid_campaign_no_freq.py
│   │   ├── validate_qed_dft_campaign.py              # Validation script
│   │   ├── ortho.xyz / meta.xyz / para.xyz           # Pristine starting geometries
│   │   ├── intermediate_scans.csv                    # (θ,φ) grid definition
│   │   ├── ortho_camp.out / para_camp.out            # Campaign log files
│   │   ├── MD_test.out                              # MD test output (empty)
│   │   ├── trajectory_sampling_direction_A.dat      # Trajectory sampling data
│   │   ├── trajectory_campaign.py                    # Trajectory optimization driver
│   │   ├── trajectory_campaign_no_freq.py           # Trajectory opt (no frequency) driver
│   │   ├── timer.dat                                # Timing data
│   │   └── progress                                 # Campaign progress log
│   ├── validate_qed_dft_campaign.py                  # Validation script (copy)
│   ├── validate_readme.md                           # Validation usage notes
│   ├── relaxed_structures.tar                       # Archive of all optimized.xyz files
│   ├── run_grid_no_freq.zip                         # Archive of runs_grid_no_freq/
│   ├── grid_scans_optimized_geoms.zip               # Archive of grid scan outputs
│   ├── timer.dat                                    # Timing data (copy)
│   └── runs_grid_no_freq/                           # 573 per-cell output directories
│       ├── grid_campaign_no_freq_{ortho,meta,para}_status.json  # Per-cell rollup (JSON)
│       ├── projected_gradient_ring_raw.csv           # Ring diagnostic data
│       ├── projected_gradient_ring_para_meta_pairwise.csv       # Pairwise ring comparison
│       ├── projected_gradient_ring_para_meta_diagnostic.png     # Ring diagnostic plot
│       ├── th63_ring_diagnostic/                     # θ=63° ring analysis (ortho↔meta)
│       ├── th72_ring_diagnostic/                     # θ=72° ring analysis
│       ├── th81_ring_diagnostic/                     # θ=81° ring analysis
│       ├── analyze_projected_gradient_ring.py        # Ring analysis script
│       ├── README_projected_gradient_ring.md        # Ring diagnostic documentation
│       └── <isomer>_th<theta>_ph<phi>/               # 573 individual cell folders
│           ├── cell.json                             # Cell manifest (isomer, θ, φ, λ-vector)
│           ├── opt_status.json                       # Optimization convergence info
│           ├── optimized.xyz                         # Final converged geometry
│           ├── opt_traj.xyz                          # Optimization trajectory
│           ├── <cell_id>_qed_ccsd_input.json         # QED-CCSD input for this cell
│           └── DONE                                  # Completion marker
│
└── qed_ccsd_on_optimized_geometries/               # ── QED-CCSD single-point campaign ──
    └── runs_grid_no_freq/                            # 573 per-cell QED-CCSD run directories
        ├── build_relaxed_qed_ccsd_scans.py          # Script: parse .out → CSV
        ├── build_grid_no_freq_energy_maps.py        # Script: merge energies + plot maps
        ├── grid_campaign_no_freq_{ortho,meta,para}_status.{csv,json}  # Rolled-up status
        ├── grid_campaign_no_freq_opt_energies.csv    # Merged opt energies (3 isomers)
        ├── ortho_meta_diff_grid_campaign_no_freq.png  # ΔE map (ortho − meta)
        ├── para_meta_diff_grid_campaign_no_freq.png   # ΔE map (para − meta)
        └── <isomer>_th<theta>_ph<phi>/               # 573 individual cell folders
            ├── cell.json                             # Cell manifest (same as relax campaign)
            ├── opt_status.json                       # Optimization info (same as relax campaign)
            ├── optimized.xyz                         # Relaxed geometry (same as relax campaign)
            ├── opt_traj.xyz                          # Optimization trajectory (same as relax campaign)
            ├── <cell_id>_qed_ccsd_input.json         # QED-CCSD input JSON
            ├── <cell_id>_qed_ccsd_input.out         # ExaChem QED-CCSD output (573 files)
            ├── validation_recompute.out              # Spot-check recompute (40 files)
            └── DONE                                  # Completion marker
```

**Note:** The `QED_DFT_Theta_Phi_Relax_Campaign/` directory also contains `.py` files at its top level (`campaign.py`, `cavity_common.py`, `grid_campaign_no_freq.py`, `ccsd_input.py`, and `validate_qed_dft_campaign.py`) that are identical copies of the files in `campaign_folder/`. These are kept for convenience as the active working directory of the campaign.

**Workflow:** The two campaigns are sequential — first QED-DFT geometries are optimized (`QED_DFT_Theta_Phi_Relax_Campaign/`), then QED-CCSD single-point energies are computed on those relaxed geometries (`qed_ccsd_on_optimized_geometries/`).

---

## The (θ, φ) Grid

The cavity field is oriented along a unit vector parameterized by spherical polar angles (θ, φ), where θ is the polar angle from the +z axis (0°–90°) and φ is the azimuthal angle in the x-y plane (0°–360°).

| Parameter | Values | Count |
|-----------|--------|-------|
| θ (polar) | 0°, 9°, 18°, 27°, 36°, 45°, 54°, 63°, 72°, 81°, 90° | 11 |
| φ (azimuthal) | 0°, 18°, 36°, 54°, 72°, 90°, 108°, 126°, 144°, 162°, 180°, ..., 342° | 20 (at θ=9–81) |

**Total unique grid points: 191**

This is not simply 11 × 20 = 220 because of pole degeneracies:

- **θ = 0° (north pole):** The field vector is (0, 0, 1) regardless of φ — only **1 point** is computed.
- **θ = 90° (equator):** By the inversion symmetry E(θ, φ) = E(180°−θ, (φ+180°) mod 360°), antipodal points are equivalent, so only φ = 0°–180° (10 points) are computed; the 180°–342° half is reconstructed by symmetry.
- **θ = 9°–81°:** All 20 φ values (0°–342°, step 18°) are computed independently.

**Grid coverage:** 1 + (9 × 20) + 10 = 191 unique points

**Symmetry expansion:** `read_csv_transform_to_full_theta_phi_grid.py` uses the inversion symmetry E(θ, φ) = E(180°−θ, (φ+180°) mod 360°) and azimuthal wrap symmetry (φ=360° → φ=0°) to reconstruct the full 21 × 20 = 420-point grid (θ: 0°–180°, φ: 0°–360°), plus the θ=0° and θ=180° polar degeneracies. This produces 440 rows in `Relaxed_QED_CCSD_Combined_Results.txt`.

---

## Computational Methods

| Parameter | Value |
|-----------|-------|
| **Molecule** | Wheland intermediates of bromo-nitrobenzene (C₆H₅Br–NO₂ derivatives), 15 atoms |
| **Isomers** | ortho, meta, para (Br-substituent position) |
| **Charge / Multiplicity** | +1 / 1 (singlet) |
| **QED-DFT optimizer** | Psi4 + `cqed_scf` ASE calculator |
| **DFT functional** | ωB97X-D (wb97x) |
| **Basis set** | 6-311G* |
| **Cavity frequency ω** | 0.06615 rad⁻¹ (≈ 929 nm) |
| **Cavity coupling λ** | 0.1 (fixed for all cells) |
| **Convergence threshold** | Projected gradient norm \|g_proj\| < 2×10⁻³ (relaxed) / 2.5×10⁻³ (full) |
| **QED-CCSD engine** | ExaChem (GPU-accelerated) |
| **Unit conversion** | 1 Hartree = 627.509 kcal/mol |

**Warm-start strategy:** Within each θ ring, φ=0° is optimized from the pristine isomer geometry; each subsequent φ is warm-started from the previous φ's converged geometry. The 11 θ rings (and 3 isomers) are independent and were run concurrently.

---

## Per-Cell Run Directory Layout

Each of the 573 grid points has a directory named `<isomer>_th<theta>_ph<phi>/` (e.g., `meta_th63_ph126/`).

### QED-DFT relaxation campaign

| File | Description |
|------|-------------|
| `cell.json` | Cell manifest: isomer, θ, φ, magnitude (0.1), λ-vector, ring/chain index |
| `opt_status.json` | Convergence info: `converged`, `attempts`, `gnorm_history`, `final_gnorm`, `final_energy_hartree`, `conv_threshold`, `wall_seconds`, `lambda_vector` |
| `optimized.xyz` | Final converged geometry (XYZ format, 15 atoms) |
| `opt_traj.xyz` | Full optimization trajectory (all BFGS restart attempts, annotated with \|g_proj\| and energy) |
| `<cell_id>_qed_ccsd_input.json` | Pre-built QED-CCSD input JSON (covers all 3 isomers at this θ, φ) |
| `DONE` | Empty marker file indicating the cell completed |

### QED-CCSD single-point campaign

| File | Description |
|------|-------------|
| `cell.json` | Same as relaxation campaign |
| `opt_status.json` | Same as relaxation campaign |
| `optimized.xyz` | Relaxed geometry (same as relaxation campaign) |
| `opt_traj.xyz` | Optimization trajectory (same as relaxation campaign) |
| `<cell_id>_qed_ccsd_input.json` | QED-CCSD input JSON |
| `<cell_id>_qed_ccsd_input.out` | **ExaChem output file** — contains SCF energy (`** Total SCF energy = ...`) and CCSD total energy (`CCSD total energy / hartree = ...`) |
| `validation_recompute.out` | Spot-check recomputation output (present in 40 cells) |
| `DONE` | Empty marker file |

---

## Data File Format Reference

### `relaxed_qed_ccsd_intermediate_scans.csv` (191 rows + header)

Primary data table. Each row = one unique (θ, φ) grid point.

| Column | Description |
|--------|-------------|
| `theta` | Polar angle in degrees (0–90) |
| `phi` | Azimuthal angle in degrees (0–342) |
| `E_SCF_ortho_int` | SCF energy (Hartree) for the ortho intermediate |
| `E_CCSD_ortho_int` | CCSD total energy (Hartree) for the ortho intermediate |
| `E_SCF_meta_int` | SCF energy (Hartree) for the meta intermediate |
| `E_CCSD_meta_int` | CCSD total energy (Hartree) for the meta intermediate |
| `E_SCF_para_int` | SCF energy (Hartree) for the para intermediate |
| `E_CCSD_para_int` | CCSD total energy (Hartree) for the para intermediate |

### `Relaxed_QED_CCSD_Combined_Results.txt` (440 rows + 2 header lines)

Symmetry-expanded grid with field vector components and CCSD energies.

| Column | Description |
|--------|-------------|
| `theta` | Polar angle in degrees (0–180, expanded) |
| `phi` | Azimuthal angle in degrees (0–360, expanded) |
| `Ex` | Field vector x-component (unit vector) |
| `Ey` | Field vector y-component (unit vector) |
| `Ez` | Field vector z-component (unit vector) |
| `Para_E` | CCSD total energy (Hartree) for the para intermediate |
| `Ortho_E` | CCSD total energy (Hartree) for the ortho intermediate |
| `Meta_E` | CCSD total energy (Hartree) for the meta intermediate |

### Status CSV (`grid_campaign_no_freq_{isomer}_status.csv`, 191 rows + header)

| Column | Description |
|--------|-------------|
| `cell_id` | e.g. `meta_th63_ph126` |
| `group_id` | Theta-ring group, e.g. `meta_th63` |
| `isomer` | `ortho` / `meta` / `para` |
| `theta` | Polar angle (degrees) |
| `phi` | Azimuthal angle (degrees) |
| `chain_index` | Position within the φ warm-start chain |
| `opt_converged` | `True`/`False` |
| `attempts` | Number of BFGS restart cycles |
| `final_gnorm` | Final projected gradient norm (Hartree/Bohr) |
| `E_opt_hartree` | Final optimized energy (Hartree) |
| `conv_threshold` | Convergence threshold used (2e-3 or 6e-3) |
| `qed_ccsd_input_written` | `True` if CCSD input JSON was generated |

### Status JSON (`grid_campaign_no_freq_{isomer}_status.json`)

Same data as the CSV, one JSON object per cell in a JSON array.

### `grid_campaign_no_freq_opt_energies.csv` (192 rows + header)

Merged QED-DFT optimization energies for all three isomers at each grid point.

| Column | Description |
|--------|-------------|
| `theta` | Polar angle (degrees) |
| `phi` | Azimuthal angle (degrees) |
| `E_opt_ortho` | QED-DFT optimized energy, ortho (Hartree) |
| `E_opt_meta` | QED-DFT optimized energy, meta (Hartree) |
| `E_opt_para` | QED-DFT optimized energy, para (Hartree) |

---

## Reproducibility & Plotting

### Regenerate `Relaxed_QED_CCSD_Combined_Results.txt` from the CSV

```bash
python read_csv_transform_to_full_theta_phi_grid.py
```

This reads `relaxed_qed_ccsd_intermediate_scans.csv` and writes:
- `Relaxed_QED_CCSD_Combined_Results.txt` (CCSD energies, symmetry-expanded)
- `Relaxed_QED_SCF_Combined_Results.txt` (SCF energies, symmetry-expanded)

### Regenerate Figure 2 plots

```bash
python plot_single_panel_ortho_para_meta.py
```

Generates:
- `relaxed_ortho_meta_diff_QED_CCSD_22.png` — ΔE (Ortho − Meta) in kcal/mol
- `relaxed_para_meta_diff_QED_CCSD_22.png` — ΔE (Para − Meta) in kcal/mol

### Regenerate QED-CCSD scan CSV from raw outputs

```bash
cd qed_ccsd_on_optimized_geometries/runs_grid_no_freq
python build_relaxed_qed_ccsd_scans.py \
    --runs-dir . \
    --grid ../../QED_DFT_Theta_Phi_Relax_Campaign/unrelaxed_qed_ccsd_intermediate_scans.csv
```

### Rebuild merged energy maps and plots

```bash
cd qed_ccsd_on_optimized_geometries/runs_grid_no_freq
python build_grid_no_freq_energy_maps.py --data-dir .
```

### QED-DFT geometries (archived)

- `relaxed_structures.tar` — contains all 573 `optimized.xyz` files from the relaxation campaign
- `grid_scans_optimized_geoms.zip` — zip archive of the grid scan outputs
- `run_grid_no_freq.zip` — zip archive of the `runs_grid_no_freq/` directory

---

## Validation & Quality Control

**All 191 grid points converged successfully** (all `opt_converged = True` in the status CSVs).

### Validation script

The QED-DFT campaign includes a validation tool:

```bash
cd QED_DFT_Theta_Phi_Relax_Campaign
python validate_qed_dft_campaign.py runs_grid_no_freq \
    --recompute-count 1 --seed 123       # Spot-check 1 cell with recompute
python validate_qed_dft_campaign.py runs_grid_no_freq --skip-incomplete   # Fast check only
```

The validator checks:
- Expected file presence (`cell.json`, `opt_status.json`, `optimized.xyz`, `*_qed_ccsd_input.json`, `DONE`)
- Consistency between `cell.json` and `opt_status.json` (matching cell IDs, λ-vector)
- Geometry consistency (XYZ atom count = 15)
- Energy consistency between JSON and XYZ comment lines

### QED-CCSD recomputation

40 cells include a `validation_recompute.out` file — spot-check recomputations at a subset of grid points to verify reproducibility of the QED-CCSD energies.

### Ring diagnostics

Three θ rings (θ = 63°, 72°, 81°) were analyzed for projected-gradient norm behavior across the φ sweep:

```
QED_DFT_Theta_Phi_Relax_Campaign/runs_grid_no_freq/th63_ring_diagnostic/
QED_DFT_Theta_Phi_Relax_Campaign/runs_grid_no_freq/th72_ring_diagnostic/
QED_DFT_Theta_Phi_Relax_Campaign/runs_grid_no_freq/th81_ring_diagnostic/
```

Each contains:
- `projected_gradient_ring_raw.csv` — per-cell gradient norms and energies
- `projected_gradient_ring_ortho_meta_pairwise.csv` — ortho↔meta comparisons at each φ
- `projected_gradient_ring_ortho_meta_diagnostic.png` — energy/gradient plots

---

## Citation

> *Cite the manuscript (DOI to be added upon publication). This dataset is hosted on Zenodo as supporting data for the associated publication.*

**Keywords:** Quantum Electrodynamics (QED) chemistry · Density Functional Theory (DFT) · Coupled Cluster Singles and Doubles (CCSD) · Cavity Quantum Electrodynamics · Wheland intermediate · geometry optimization · energy scan · bromo-nitrobenzene
