# qed_dft_relaxed - Relaxed QED-DFT energy tables

This folder contains per-lambda relative-energy tables from geometries relaxed
under QED-DFT while holding the cavity orientation fixed at the target
direction. Each row is one coupling magnitude.

## Files

| File | Pair | Direction | Rows |
| --- | --- | --- | --- |
| `ortho-meta_70_31_hartree.csv` | ortho - meta | `(theta=70, phi=31)` | 5 |
| `para-meta_65_78_hartree.csv` | para - meta | `(theta=65, phi=78)` | 5 |
| `availability.md` | both pairs | both directions | Human-readable availability summary. |

## Columns

The CSVs use the same leading metadata columns as the other method summaries.

| Column | Meaning |
| --- | --- |
| `theta`, `phi` | cavity direction in degrees |
| `Ex`, `Ey`, `Ez` | Cartesian components of the lambda vector |
| `lambda_magnitude` | `|lambda|` in a.u. |
| `E_<isomer>_Hartrees` | absolute electronic energy in Hartree |
| `zpe_<isomer>_Hartrees` | harmonic zero-point energy in Hartree |
| `dE_<pair>_raw_kcal/mol` | electronic relative energy, `E(A) - E(B)`, in kcal/mol |
| `dE_<pair>_zpe_kcal/mol` | ZPE-corrected relative energy, `(E_A + ZPE_A) - (E_B + ZPE_B)`, in kcal/mol |
| `<isomer>_converged` | QED-DFT optimization convergence flag |
| `<isomer>_has_freq` | whether the frequency/ZPE calculation is present |

## Method notes

- Geometry type: relaxed at each `(isomer, direction, |lambda|)` point.
- Relative energies are already converted to kcal/mol using 627.509
  kcal/mol per Hartree.
- These are QED-DFT energies, not QED-CCSD energies. The matching QED-CCSD
  relaxed-geometry summaries are in `../QED_CCSD/summary/`.
