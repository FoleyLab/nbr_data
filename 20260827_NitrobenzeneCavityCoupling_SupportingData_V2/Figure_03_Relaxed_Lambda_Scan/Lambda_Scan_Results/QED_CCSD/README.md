# QED_CCSD - QED-CCSD and QED-HF calculations

This folder contains QED-CCSD(2,2) and QED-HF calculations performed on both
QED-DFT-relaxed and unrelaxed geometries. The long calculation output files are
kept unchanged in `outputs/`; reduced energy tables are in `summary/`.

## Layout

| Folder | Contents |
| --- | --- |
| `outputs/` | Original long `.out` files from individual QED-CCSD jobs, grouped by geometry and cavity direction. |
| `summary/` | Small reduced tables for QED-CCSD(2,2), QED-HF, and duplicate pQED summaries used for comparisons. |

## Directions and pairs

| Direction | Primary pair | Geometry sets |
| --- | --- | --- |
| `(theta=70, phi=31)` | `ortho - meta` | relaxed and unrelaxed |
| `(theta=65, phi=78)` | `para - meta` | relaxed and unrelaxed |

The relaxed calculations use geometries optimized with QED-DFT at the same
cavity direction and coupling magnitude. The unrelaxed calculations use fixed
reference geometries.

## Summary conventions

New `.csv` summaries are the preferred machine-readable files. They preserve
the energies from the original `.txt` summaries and add uniform metadata:

| Column | Meaning |
| --- | --- |
| `method` | `QED-CCSD(2,2)` or `QED-HF` |
| `geometry` | `relaxed` or `unrelaxed` |
| `theta`, `phi` | cavity direction in degrees |
| `Ex`, `Ey`, `Ez` | Cartesian components of the lambda vector |
| `lambda_magnitude` | `|lambda|` in a.u. |
| `E_<isomer>_Hartrees` | absolute electronic energies in Hartree |
| `dE_<pair>_kcal/mol` | relative energy, `E(A) - E(B)`, in kcal/mol |

The original `.txt` summaries are left in place for provenance.
