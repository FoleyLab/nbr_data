# Lambda_Scan_Results - Cavity-coupled isomer energetics

This folder collects relative energetics for bromonitrobenzene intermediates
under strong directional cavity coupling. The data come from heterogeneous
ab initio cavity-QED methods: QED-DFT, QED-CCSD(2,2), QED-HF, and pQED.

The main comparisons are:

| Pair | Direction | Primary relative energy |
| --- | --- | --- |
| ortho/meta | `(theta=70, phi=31)` | `ortho - meta` |
| para/meta | `(theta=65, phi=78)` | `para - meta` |

## Layout

| Folder | Contents |
| --- | --- |
| `qed_dft_relaxed/` | QED-DFT energies at geometries relaxed under the target cavity coupling. |
| `qed_dft_unrelaxed/` | QED-DFT energies at fixed unrelaxed geometries. |
| `QED_CCSD/` | QED-CCSD(2,2), QED-HF, and original long calculation outputs. |
| `pqed/` | pQED unrelaxed scans for 49 electrons, `Nph=3/10`, with and without coherent-state transform. |
| `scripts/` | Plotting script and generated comparison figures. |
| `docs/` | Presentation material derived from these summaries. |

Each subfolder has its own `README.md` describing the files it contains.

## Data conventions

Summary CSVs use consistent leading metadata columns where possible:

| Column | Meaning |
| --- | --- |
| `method` | Present when a file may be mixed or generated from method-specific text summaries. |
| `geometry` | `relaxed` or `unrelaxed`, where applicable. |
| `theta`, `phi` | cavity direction in degrees. |
| `Ex`, `Ey`, `Ez` | Cartesian components of the lambda vector. |
| `lambda_magnitude` | `|lambda|` in a.u. |
| `E_<isomer>_Hartrees` | absolute electronic energy in Hartree. |
| `zpe_<isomer>_Hartrees` | zero-point energy in Hartree, where available. |
| `dE_<pair>_..._kcal/mol` | relative energy in kcal/mol. |

Relative energies use the sign convention `E(A) - E(B)`. The two target
columns are `dE_ortho_meta*_kcal/mol` and `dE_para_meta*_kcal/mol`.

## Geometry and method notes

- QED-DFT relaxed data use geometries optimized at each target direction and
  coupling magnitude.
- QED-DFT unrelaxed, QED-CCSD unrelaxed, and pQED use fixed unrelaxed
  geometries.
- QED-CCSD relaxed calculations reuse the QED-DFT-relaxed geometries.
- pQED is only available for unrelaxed geometries.
- The long `.out` calculation files in `QED_CCSD/outputs/` are left unchanged
  and documented rather than reorganized or parsed in place.

## Plotting

Use `scripts/plot_isomer_energies.py` to regenerate comparison plots:

```bash
cd scripts
python plot_isomer_energies.py
```

The default script configuration generates comparisons for QED-CCSD relaxed
vs. unrelaxed, QED-DFT relaxed vs. unrelaxed, pQED vs. QED-CCSD, and pQED vs.
QED-DFT for the target isomer pairs and directions.
