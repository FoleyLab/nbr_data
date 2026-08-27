# qed_dft_unrelaxed — Unrelaxed QED-DFT scans

Relative-energy scans of the two isomer pairs vs. cavity coupling |λ| using
**unrelaxed (fixed gas-phase) geometries**. One row per λ (0.02–0.10 a.u.).
Electronic energies are in Hartree; the difference is reported in kcal/mol.
Columns follow the pQED scan-file convention.

## Files

| File | Pair | Direction | Rows |
| --- | --- | --- | --- |
| `ortho_meta_dir70_31_scan_qed_dft_no_relax.csv` | ortho − meta | (θ=70°, φ=31°) | 5 |
| `para_meta_dir65_78_scan_qed_dft_no_relax.csv` | para − meta | (θ=65°, φ=78°) | 5 |

## Columns

| Column | Meaning |
| --- | --- |
| `theta`, `phi` | cavity direction (degrees) |
| `Ex`, `Ey`, `Ez` | Cartesian components of the λ vector (a.u.) |
| `lambda_magnitude` | \|λ\| (a.u.), the norm of (Ex, Ey, Ez) |
| `E_ortho_Hartrees` / `E_para_Hartrees` | electronic energy of the named isomer (Hartree) |
| `E_meta_Hartrees` | electronic energy of meta (Hartree) |
| `dE_ortho_meta_kcal/mol` / `dE_para_meta_kcal/mol` | energy difference, A − B, in kcal/mol |

## Method notes

- Geometry type: **unrelaxed** — the same gas-phase geometry is used at every
  λ for a given isomer/direction.
- λ components are generated from the spherical direction:
  `Ex = λ·sinθ·cosφ`, `Ey = λ·sinθ·sinφ`, `Ez = λ·cosθ`.
- No ZPE correction: not available/meaningful for unrelaxed geometries.
- Relative energies are already in kcal/mol (× 627.509); do not convert the
  `dE_*` columns again.
