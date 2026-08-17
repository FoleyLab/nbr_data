
**Data for:** Figure S1 — [Target Calculation Details]

**System:** Three tautomeric Wheland intermediates of bromo-nitrobenzene — **ortho**, **meta**, and **para** (relative to the Br substituent) — each (15 atoms, charge +1, doublet-derived singlet) evaluated across a spherical polar (θ, φ) grid. Multiple levels of theory are included including QED-DFT and QED-CCSD.

---

## Table of Contents
1. [Quick Start](#quick-start)
2. [Directory Overview](#directory-overview)
3. [The (θ, φ) Grid](#the-grid)
4. [Computational Methods](#computational-methods)
5. [Data File Format Reference](#data_file_format_reference)
6. [Reproducibility & Plotting](#reproducibility--plotting)

---

## Quick Start

For consolidated data and results:

| File | Description |
|------|-------------|
| `unrelaxed_QED_CCSD_theta_phi_scan.txt` | Primary QED-CCSD scan data |
| `unrelaxed_pQED_theta_phi_scan.txt` | pQED dataset (standard) |
| `unrelaxed_pQED_CS_theta_phi_scan.txt` | pQED dataset with CS correction |
| `unrelaxed_qeddft_theta_phi_scan.txt` | QED-DFT data |
| `[suffix]_ortho_meta_*.png` | Energy difference maps for each respective dataset |

To generate the plots using the script:

```bash
python plot_single_panel_ortho_para_meta.py
Directory Overview
Figure_S1_Unrelaxed_Energy_Scans/
├── README.md                                          # This file
├── plot_single_panel_ortho_para_meta.py             # Script: generate difference plots for all 4 datasets
├── unrelaxed_QED_CCSD_theta_phi_scan.txt            # Data: QED-CCSD results
├── unrelaxed_pQED_theta_phi_scan.txt                # Data: pQED results
├── unrelaxed_pQED_CS_theta_phi_scan.txt             # Data: pQED + CS result
└── unrelaxed_qeddft_theta_phi_scan.txt              # Data: QED-DFT results
The (θ, φ) Grid
The cavity field is oriented along a unit vector parameterized by spherical polar angles (θ, φ), where θ is the polar angle from the +z axis (0°–90°) and φ is the azimuthal angle in the x-y plane (0°–360°).
Parameter	Values
θ (polar)	0°, 9°, 18°, 27°, 36°, 45°, 54°, 63°, 72°, 81°, 90°
φ (azimuthal)	0°, 18°, 36°, ..., 342°
Total unique grid points: 191
Computational Methods
Parameter
Molecule
Isomers
Charge / Multiplicity
QED-DFT Theory
QED-CCSD Engine
Unit conversion
Data File Format Reference
Data Input Files (*.txt)
These files contain the raw calculation data for each point on the grid:
- theta (col 0)
- phi (col 1)
- e_para (col 5)
- e_ortho (col 6)
- e_meta (col 7)
Reproducibility & Plotting
To generate the publication-quality difference plots for all four dataset variants:
python plot_single_panel_ortho_para_meta.py
This will output 8 total PNG files, grouped in pairs by data source:
1. QED_CCSD (Standard)
2. pQED 
3. pQED_CS
4. qeddf (DFT Reference)
