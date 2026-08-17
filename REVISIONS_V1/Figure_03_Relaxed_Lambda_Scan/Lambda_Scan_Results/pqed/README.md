# pqed — Parameterized QED (pQED) scans

Unrelaxed-geometry parameterized QED scans of all three bromonitrobenzene
isomers at two cavity directions. Each file covers λ = 0.02–0.10 a.u. (5 rows)
and reports energies of all isomers (Hartree) plus differences for every
isomer pair (kcal/mol).

## Files

8 files = 2 directions × 2 photon-state counts × 2 coherent-state options.

| Direction | Nph | Non-CS | Coherent-state (CS) |
| --- | --- | --- | --- |
| (70°, 31°) | 3 | `pqed_49_3_dir70_31_scan.csv` | `pqed_49_3_dir70_31_scan_CS.csv` |
| (70°, 31°) | 10 | `pqed_49_10_dir70_31_scan.csv` | `pqed_49_10_dir70_31_scan_CS.csv` |
| (65°, 78°) | 3 | `pqed_49_3_dir65_78_scan.csv` | `pqed_49_3_dir65_78_scan_CS.csv` |
| (65°, 78°) | 10 | `pqed_49_10_dir65_78_scan.csv` | `pqed_49_10_dir65_78_scan_CS.csv` |

All use **Nel = 49** electrons.

## Columns

| Column | Meaning |
| --- | --- |
| `theta`, `phi` | cavity direction (degrees) |
| `Ex`, `Ey`, `Ez` | Cartesian components of the λ vector (a.u.) |
| `lambda_magnitude` | \|λ\| (a.u.) |
| `E_ortho_Hartrees`, `E_meta_Hartrees`, `E_para_Hartrees` | pQED energies (Hartree) |
| `dE_ortho_meta_kcal/mol`, `dE_para_meta_kcal/mol`, `dE_ortho_para_kcal/mol` | pair differences (kcal/mol) |

## Method notes

- pQED: parameterized QED with 49 electrons and a truncated photonic basis of
  `Nph` states; the `_CS` variant applies a coherent-state (polaron) transform
  of the photon mode before truncation.
- Energies (~ −3007.9 Ha) are the unrelaxed-geometry pQED approximation and
  differ from the QED-CCSD(2,2) values in `../../../QED_CCSD/summary/`.
- kcal/mol conversion factor used: 627.509.

## Duplication warning

These CSVs duplicate `../QED_CCSD/summary/unrelaxed_dir_*_pqed_49_*.csv`.
The copies in this `pqed/` folder are the clearer method-level copies used by
`../scripts/plot_isomer_energies.py`; keep the two locations in sync if either
copy is regenerated.
