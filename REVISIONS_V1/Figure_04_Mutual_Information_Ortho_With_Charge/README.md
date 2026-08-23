# Mutual-Information & Resonance-Reweighting Analysis (QED-CI / MPS)

**Data for:** Figure 4 — Mutual Information of the Ortho Wheland Intermediate, with Charge Analysis

**System:** Bromo-nitrobenzene Wheland intermediates (**ortho**, **meta**, **para**) and the unsubstituted nitrobenzene reference (**unsubs**), each examined under the two cavity-field orientations used throughout this work — **A** (θ = 70°, φ = 31°) and **D** (θ = 65°, φ = 78°) — plus a **no_coupling** (field-free) baseline. Pairwise mutual information among the pz orbitals of the π system quantifies how the cavity field reweights the resonance/conjugation pattern.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Directory Overview](#directory-overview)
3. [Computational Methods](#computational-methods)
4. [Data File Format Reference](#data-file-format-reference)
5. [Reproducibility & Plotting](#reproducibility--plotting)
6. [Citation](#citation)

---

## Quick Start

Final figures live at the top level of this folder:

| File | Description |
|------|-------------|
| `Mutual_Information.png` | Mutual-information summary figure |
| `cavity_resonance_reweighting_with_charge.png` | Cavity-induced resonance reweighting, with charge analysis |
| `Resonance_Structures.png` | Resonance-structure scheme |
| `Ortho.pdf` / `Plot_Ortho_A.pdf` / `Plot_Ortho_diff_A.pdf` | Ortho-isomer mutual-information plots (orientation A and A-vs-reference difference) |

Per-isomer difference plots (`Plot_<Isomer>_diff_<A|D>.pdf`) are co-located with their input data under `Mut_info_plots/`.

---

## Directory Overview

```
Figure_04_Mutual_Information_Ortho_With_Charge/
├── README.md                                   # This file
├── Mutual_Information.png                      # Final figure
├── cavity_resonance_reweighting_with_charge.png# Final figure
├── Resonance_Structures.png                    # Final figure (scheme)
├── Ortho.pdf / Plot_Ortho_A.pdf / Plot_Ortho_diff_A.pdf   # Final ortho plots
│
├── Mut_info_calculations/NitroBenzene/         # ── Raw QED-CI/MPS run directories ──
│   ├── ortho/  meta/  para/  unsubs/           # One folder per species
│   │   ├── A_orientation/                      # Cavity field along direction A
│   │   ├── D_orientation/                      # Cavity field along direction D
│   │   └── no_coupling/                        # Field-free baseline
│   │       ├── FCIDUMP_el                      # Electronic integrals (FCIDUMP format)
│   │       ├── FCIDUMP_int                     # Electron–photon interaction integrals
│   │       ├── FCIDUMP_ph                      # Photon terms
│   │       ├── molmps.inp                      # MPS solver input
│   │       ├── molmps.inp.out / .tmp           # MPS solver output / scratch
│   │       ├── correls_2.moltools              # Two-orbital correlation output
│   │       ├── QED-CI_test.out                 # QED-CI run log
│   │       ├── qedci                           # Cluster submission script (SLURM)
│   │       ├── psi4_<jobid>/                   # Psi4 scratch from integral generation
│   │       └── slurm-<jobid>.out               # Cluster scheduler logs
│
└── Mut_info_plots/                             # ── Analysis / plotting working dirs ──
    ├── Ortho/ Meta/ Para/                      # Br-substituted intermediates
    │   ├── A_orientation/                      # Calc artifacts (as above) +
    │   │   ├── xyz                             # Geometry + pz-orbital index per atom
    │   │   ├── mut_info                        # N×N pairwise mutual-information matrix
    │   │   ├── script.py                       # Plotting script
    │   │   └── Plot_<Isomer>_diff_A.pdf        # Difference plot (vs. reference)
    │   └── D_orientation/                      # Same layout for direction D
    └── Unsub/                                  # Unsubstituted nitrobenzene
        ├── A_orientation/  D_orientation/      # Per-orientation working dirs
        ├── mut_info / xyz / n_alpha            # Reference data (no charge overlay)
        ├── script.py / script_Nocharge.py      # Plotting scripts
        └── Plot_Unsubst.pdf / Plot_Unsubst_charges.pdf
```

**Note:** The `Mut_info_plots/<Isomer>/<orientation>/` directories each contain a full copy of the calculation inputs (FCIDUMP files, `molmps.*`, etc.) alongside the extracted `mut_info`/`xyz` data and the plotting script — they are self-contained working directories from the analysis machine. The `psi4_*` directories and `slurm-*.out` files are compute-cluster artifacts retained for provenance.

---

## Computational Methods

| Parameter | Value |
|-----------|-------|
| **Species** | ortho / meta / para Wheland intermediates + unsubstituted nitrobenzene |
| **Integral generation** | Psi4 (FCIDUMP format) |
| **Correlated solver** | QED-CI / MPS (`molmps` input; `qed-ci` code, run on a CPU cluster via SLURM) |
| **Cavity orientations** | Direction A (θ = 70°, φ = 31°), direction D (θ = 65°, φ = 78°), and field-free (`no_coupling`) |
| **Analysis quantity** | Pairwise mutual information among pz orbitals of the π system (`correls_2.moltools` → `mut_info`) |

---

## Data File Format Reference

### `mut_info`

Symmetric N×N whitespace-delimited matrix of pairwise mutual information between pz orbitals (N = 10 for the Br-substituted intermediates: 6 ring C + Br + N + 2 O). Row/column ordering follows the pz-orbital indices given in the companion `xyz` file.

### `xyz`

Whitespace-delimited, no header: `element  x  y  z  pz_index`, where `pz_index` (4th column) labels the pz orbital carried by that atom, i.e. the row/column index used in `mut_info`.

### FCIDUMP files

| File | Contents |
|------|----------|
| `FCIDUMP_el` | Standard electronic one- and two-electron integrals |
| `FCIDUMP_int` | Electron–photon (bilinear coupling) integrals |
| `FCIDUMP_ph` | Photon contribution terms |

---

## Reproducibility & Plotting

The mutual-information plots are generated per isomer/orientation from the working directories:

```bash
cd Mut_info_plots/Ortho/A_orientation
python script.py xyz mut_info
```

`script.py` takes two positional arguments (geometry file, mutual-information file) and also reads the same filenames from the parent directory (`../`) to build the comparison/difference panels — run it from within the orientation directory as shown. The unsubstituted reference uses `script.py` / `script_Nocharge.py` in `Mut_info_plots/Unsub/`.

---

## Citation

> *Cite the manuscript (DOI to be added upon publication). This dataset is hosted on Zenodo as supporting data for the associated publication.*

**Keywords:** Quantum Electrodynamics (QED) chemistry · QED-CI · MPS/DMRG · mutual information · resonance structures · Wheland intermediate · bromo-nitrobenzene
