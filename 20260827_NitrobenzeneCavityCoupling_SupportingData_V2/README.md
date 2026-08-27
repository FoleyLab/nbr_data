# Supporting Data — Cavity-QED Control of Bromo-Nitrobenzene Wheland Intermediates (Revision V1)

This repository contains the computational data, analysis scripts, and final figures supporting the revised manuscript on tuning the relative stability of bromo-nitrobenzene **Wheland intermediates** (ortho, meta, para) with a strongly coupled optical cavity. Cavity-field **orientation** (spherical polar angles θ, φ) and **coupling strength** (λ) are scanned and the isomer energetics are evaluated at multiple levels of theory (QED-DFT, QED-CCSD, pQED/Pauli–Fierz), complemented by ab initio molecular dynamics and mutual-information analysis.

---

## Repository Map

Each folder is self-contained: it has its own `README.md` with a Quick Start, the primary data tables, and the scripts needed to regenerate its figure(s).

| Folder | Manuscript item | Contents |
|--------|-----------------|----------|
| `Figure_01_Scheme/` | Figure 1 | Scheme of nitrobenzene and the ortho/meta/para intermediates on a common coordinate system (`Scheme_1.png`) |
| `Figure_02_Relaxed_Theta_Phi_Scan/` | **Figure 2** | Central result: QED-DFT geometry relaxations + QED-CCSD single points on a 191-point (θ, φ) cavity-orientation grid (573 cells); energy-difference maps |
| `Figure_03_Relaxed_Lambda_Scan/` | Figure 3 | Coupling-strength (λ) scans at fixed directions; relaxed vs. unrelaxed QED-CCSD comparisons |
| `Figure_04_Mutual_Information_Ortho_With_Charge/` | Figure 4 | QED-DMRG mutual-information and resonance-reweighting analysis (directions A, D*, and field-free) |
| `Figure_05_Molecular_Dynamics_on_Relaxed_Surface/` | Figure 5 | AIMD trajectories (directions A, D*) overlaid on the relaxed energy-difference surfaces |
| `Figure_06_Molecular_Dynamics_Relaxed_Timeseries/` | Figure 6 | ΔE timeseries and basin dwell-time analysis along the trajectories |
| `Figure_S1_Unrelaxed_Theta_Phi_Scans/` | Figure S1 | (θ, φ) scans on unrelaxed geometries at four theory levels (pQED, pQED+CS, QED-DFT, QED-CCSD) |
| `Figure_S2_Unrelaxed_Energy_Decomposition/` | Figure S2 | pQED energy decomposition (E_el/E_ph/E_blc/E_dse) from ChronusQ EOM-CCSD data |
| `Figure_S7_Snapshots/` | Figure S7 | Rendered structure snapshots (direction-A trajectory; Wheland intermediates) |
| `Figure_S8_Nitrobenzene_Energy_Fluctuations/` | Figure S8 | Cavity vs. cavity-free potential-energy fluctuations along the direction-A trajectory |
| `Table_S1_Triples_Correction/` | Table S1 | ORCA CCSD(T) and DLPNO-CCSD(T) triples-correction checks for the three isomers |

---
* Direction D is referred to as Trajectory B in the manuscript.
## Common Computational Details

Shared across folders unless a folder's README states otherwise:

| Parameter | Value |
|-----------|-------|
| **Wheland intermediates** | C₆H₅Br–NO₂ (ortho/meta/para), 15 atoms, charge +1, singlet |
| **Trajectory species** | Nitrobenzene, 14 atoms |
| **DFT functional / basis** | ωB97X-D / 6-311G* |
| **Cavity frequency ω** | 0.06615 a.u. (≈ 688 nm) |
| **Coupling λ** | 0.1 for the Figure 2 orientation scan; scanned 0–0.10 in Figures 3 / S2 |
| **Unit conversion** | 1 Hartree = 627.509 kcal/mol |
| **Orientation grid** | θ ∈ [0°, 90°] (11 values), φ ∈ [0°, 342°] (20 values per θ ring); 191 unique points, expanded to a 440-row full grid via the inversion symmetry E(θ, φ) = E(180°−θ, (φ+180°) mod 360°) — see the Figure_02 README |

### Recurring cavity-field directions

Two field orientations recur throughout the repository (in folder/file names as `dir70_31` / `dir65_78`, `direction_A` / `direction_D`, or `A_orientation` / `D_orientation`):

| Label | (θ, φ) | Associated isomer pair |
|-------|--------|------------------------|
| **Direction A** | (70°, 31°) | ortho vs. meta |
| **Direction D** | (65°, 78°) | para vs. meta |

*Note*: **Direction D** is referred to as **Trajectory B** in the manuscript.
---

## Software

| Code | Used for |
|------|----------|
| Psi4 + `cqed_scf` (ASE) | QED-DFT optimizations, AIMD, energy/gradient evaluations |
| ExaChem | GPU-accelerated QED-CCSD single points (Figure 2) |
| ChronusQ | EOM-CCSD reference data for the pQED Hamiltonians (Figure S2) |
| ORCA | CCSD(T) / DLPNO-CCSD(T) triples corrections (Table S1) |
| MolMPS | QED-DMRG Mutual-information calculations (Figure 4) |
| Python (numpy, pandas, matplotlib) | All analysis and plotting |

---

## How to Use This Repository

1. Start with the manuscript figure of interest in the [Repository Map](#repository-map) and open that folder's `README.md`.
2. Each README's **Quick Start** table points directly to the primary data files and final figures — no need to dig into raw calculation output.
3. Regeneration commands (CSV → combined grids → plots) are listed in each folder's **Reproducibility & Plotting** section.

## Validation Summary

- **Figure 2:** all 191 (θ, φ) grid points converged (all isomers); 40 cells spot-checked with QED-CCSD recomputations; per-cell validation script included.
- **Figure S2:** dipole-moment convention of the EOM-CCSD data verified against SCF/ADC2 references (`dipole_convention_check.py`).
- **Table S1:** triples corrections computed at two levels (CCSD(T) and DLPNO-CCSD(T)).

---

## Citation

> *Cite the manuscript (DOI to be added upon publication). This dataset is hosted on Zenodo as supporting data for the associated publication.*

**Keywords:** cavity quantum electrodynamics · QED-DFT · QED-CCSD · Pauli–Fierz · polariton chemistry · Wheland intermediate · bromo-nitrobenzene · molecular dynamics · mutual information
