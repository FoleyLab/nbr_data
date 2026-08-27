# Unrelaxed pQED Energy Decomposition from ChronusQ EOM-CCSD Data

**Data for:** Figure S2 — Pauli–Fierz (pQED) Energy Decomposition on Unrelaxed Geometries

**System:** Bromo-nitrobenzene Wheland intermediates (**ortho**, **meta**, **para**) on their pristine (unrelaxed) geometries. A Pauli–Fierz (PF) Hamiltonian is built from ChronusQ EOM-CCSD data and diagonalized to obtain polaritonic ground-state energies, which are then decomposed into electronic, photon, bilinear-coupling, and dipole-self-energy components as functions of coupling strength λ at fixed cavity-field directions.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Directory Overview](#directory-overview)
3. [Computational Methods](#computational-methods)
4. [Outputs & Naming Convention](#outputs--naming-convention)
5. [Validation: Dipole-Convention Check](#validation-dipole-convention-check)
6. [Citation](#citation)

---

## Quick Start

| File | Description |
|------|-------------|
| `pqed_direct_decomposition_ortho_meta_dir70_31_49_10_cs.png` | Ortho − Meta decomposition, direction (θ=70°, φ=31°), Nel=49 / Nph=10, CS Hamiltonian |
| `pqed_direct_decomposition_ortho_meta_dir70_31_49_3_cs.png` | Same pair/direction, Nph=3 |
| `pqed_direct_decomposition_para_meta_dir65_78_49_10_cs.png` | Para − Meta decomposition, direction (θ=65°, φ=78°), Nel=49 / Nph=10, CS Hamiltonian |
| `pqed_direct_decomposition_para_meta_dir65_78_49_3_cs.png` | Same pair/direction, Nph=3 |
| `pqed_direct_decomposition_para_meta_dir65_78_49_2_cs.png` | Same pair/direction, Nph=2 |

To regenerate the decomposition plots (requires the `p4dev` conda environment):

```bash
conda run -n p4dev python plot_pqed_decomp_direct.py
```

---

## Directory Overview

```
Figure_S2_Unrelaxed_Energy_Decomposition/
├── README.md                              # This file
├── pf_isomer_scan.py                      # PF Hamiltonian (θ,φ) scan driver
├── plot_pqed_decomp_direct.py             # Direct λ-scan decomposition plotter
├── plot_single_panel_ortho_para_meta.py   # Shared single-panel ΔE map plotter
├── pqed_direct_decomposition_*.png        # Output figures (see naming convention)
├── __pycache__/                           # Python bytecode cache (ignore)
└── ChronusQData/                          # ── Raw EOM-CCSD reference data ──
    ├── ortho.h5 / meta.h5 / para.h5       # ChronusQ EOM-CCSD data (HDF5)
    ├── ortho.out / meta.out / para.out    # ChronusQ output files
    ├── ortho_dipole_data.txt              # SCF + ADC2 reference dipoles (ortho)
    └── dipole_convention_check.py         # Dipole-convention diagnostic
```

---

## Computational Methods

| Parameter | Value |
|-----------|-------|
| **Reference data** | ChronusQ EOM-CCSD (`ChronusQData/*.h5`) |
| **Hamiltonian** | Pauli–Fierz, built per (θ, φ) from parsed EOM-CCSD energies and dipole matrices |
| **Photon frequency ω** | 0.066148 Hartree |
| **Coupling magnitude** | λ = 0.1 for the (θ, φ) scan; λ = 0.0–0.10 (step 0.01) for decomposition plots |
| **Basis size knobs** | `NUM_ELECTRONIC_STATES` (Nel), `NUM_FOCK_STATES` (Nph); scan: Nel=49, Nph=10 |
| **(θ, φ) scan grid** | θ ∈ [0°, 180°], φ ∈ [0°, 360°], 90 × 90 points (`pf_isomer_scan.py`) |
| **Nuclear dipoles** | Hardcoded per isomer in `pf_isomer_scan.py` |
| **Isomer labels** | Directions identified as `dir<theta>_<phi>`, e.g. `dir70_31` = (θ=70°, φ=31°) = direction A; `dir65_78` = (θ=65°, φ=78°) = direction D |

### Two Hamiltonian flavors

Both are built and diagonalized side by side for every grid point and isomer:

- **non-CS** — bilinear coupling and dipole self-energy operators built from the **total** (electronic + nuclear) dipole, so both contribute to the `d` operator.
- **CS** — the coherent-state-transformed Hamiltonian.

### `pf_isomer_scan.py` outputs

Written under `Nel_<N>_Nph_<M>/` as fixed-width text files:

| File | Contents |
|------|----------|
| `isomer_..._total_energies.dat` | Total energy per isomer, non-CS |
| `isomer_..._total_energies_CS.dat` | Total energy per isomer, CS |
| `isomer_..._energy_decomposition.dat` | E_el / E_ph / E_blc / E_dse per isomer, non-CS |
| `isomer_..._energy_decomposition_CS.dat` | Same, CS |
| `isomer_..._differences.dat` | Per-isomer (non-CS − CS) differences |

(These scan outputs are **not** stored in this folder; the `isomer_Nel_49_Nph_10_total_energies.dat` grid is referenced by the Figure_05/Figure_06 analysis scripts.)

### `plot_pqed_decomp_direct.py`

Computes the pQED component expectation values **directly from the EOM-CCSD data** (parsing `ChronusQData/*.h5` via `pf_isomer_scan.py`'s parsers) rather than reading precomputed decomposition CSVs. Key configuration (top of file): `ISOMER_A` / `ISOMER_B` (default para/meta), `THETA` / `PHI` (default 65°/78°), `NUM_ELECTRONIC_STATES` = 49, `NUM_FOCK_STATES`, `LAMBDA_MAGNITUDES` = 0.0–0.10, `TRANSFORMATION` = `"non-cs"` / `"cs"` / `"both"`, and `SUBTRACT_ZERO_BASELINE` (shifts curves to their λ=0 baseline).

---

## Outputs & Naming Convention

```
pqed_direct_decomposition_<isomerA>_<isomerB>_dir<theta>_<phi>_<Nel>_<Nph>[_cs].png
```

- `dir70_31` → ortho–meta pair at (θ=70°, φ=31°) (direction A)
- `dir65_78` → para–meta pair at (θ=65°, φ=78°) (direction D)
- Trailing `_cs` → coherent-state-transformed Hamiltonian (its absence denotes the non-CS variant when `TRANSFORMATION = "both"`)

PNG and PDF are both written when `SAVE_PNG` / `SAVE_PDF` are enabled (only the PNGs are retained in this folder).

---

## Validation: Dipole-Convention Check

`ChronusQData/dipole_convention_check.py` numerically tests which dipole-moment convention the ChronusQ EOM-CCSD `.h5` data use, comparing against SCF (ground-state) and ADC2 (excited/transition) reference data in `ortho_dipole_data.txt`. Hypotheses tested for the ground/excited **state** dipoles:

| ID | Hypothesis |
|----|------------|
| H1 | h5 state dipole = Total (electronic + nuclear), same sign |
| H2 | h5 state dipole = Electronic only, same sign |
| H3 | h5 state dipole = −Electronic only (sign-flipped) |
| H4 | h5 state dipole = Nuclear only |
| H5 | h5 state dipole − Nuclear = Electronic (i.e. h5 = Total, possibly with a different nuclear-dipole reference) |

Transition dipoles are purely electronic in both codes and are compared directly (by magnitude, since excited-state phases are arbitrary), both symmetrized bra/ket-averaged (as `pf_isomer_scan.py` uses) and raw. Run it to print the full comparison report:

```bash
python ChronusQData/dipole_convention_check.py
```

---

## Citation

> *Cite the manuscript (DOI to be added upon publication). This dataset is hosted on Zenodo as supporting data for the associated publication.*

**Keywords:** Pauli–Fierz Hamiltonian · pQED · EOM-CCSD · ChronusQ · energy decomposition · coherent-state transformation · Wheland intermediate · bromo-nitrobenzene
