# nbr_data

A data repository accompanying the publication: 
> **"Strong Light–Matter Coupling as a Photonic Substituent: Correlation-Enhanced Control of Regioselectivity in Nitrobenzene"** by Manderna et al.

---

## Repository Contents

| Directory | Description | Key Components |
| :--- | :--- | :--- |
| `COORDINATES/` | Structural coordinates for all studied species. | `.xyz` files, [Exachem](https://github.com/ExaChem/exachem) input generation script for QED-CCSD calculations. |
| `Energy_Scan_Data/` | Energy data of Wheland intermediates as a function of orientation ($\theta$, $\phi$) in a cavity field. | Formatted data, processing scripts (QED-CCSD and QED-DFT levels). |
| `Mutual_Information_Data/` | Data for mutual information analysis of Wheland intermediates (inside and outside cavity). | Raw outputs, processed data from QED-DMRG calculations. |
| `Molecular_Dynamics_Data/` | Data and scripts for the QED-DFT-driven molecular dynamics simulations. | `.xyz` trajectory outputs, processing scripts. |
| `nitrobenzene_cqed_scf_scripts/` | Driver scripts for various [QED-SCF](https://github.com/FoleyLab/cqed-scf) calculations | `.py` and `.ipynb` scripts. |

---

## Getting Started

To clone this repository, run `git clone https://github.com/FoleyLab/nbr_data.git`
