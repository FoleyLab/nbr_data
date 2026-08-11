# qed_dft_relaxed — Relaxed QED-DFT energy tables

Per-lambda relative-energy tables for the two isomer pairs, computed from
geometries **relaxed under QED-DFT** for each (isomer, cavity direction,
|λ|) combination. Column conventions follow the pQED scan files: absolute
electronic energies and ZPEs in **Hartree**, relative energy differences in
**kcal/mol** (Hartree × 627.509). Figures are produced by
`scripts/plot_isomer_energies.py`.

Produced by `../../analyze.py`, which reads `opt_status.json`
(`final_energy_hartree`) and `frequencies.json` (`zpe_hartree`) from the
per-run `runs/` folder.

## Files

| File | Source | Contents |
| --- | --- | --- |
| `ortho-meta_70_31_hartree.csv` | `../../analyze.py` | ortho − meta at (θ=70°, φ=31°), 5 rows (λ = 0.02–0.10) |
| `para-meta_65_78_hartree.csv` | `../../analyze.py` | para − meta at (θ=65°, φ=78°), 5 rows (λ = 0.02–0.10) |
| `availability.md` | `../../analyze.py` | which λ points have raw / +ZPE differences and any missing data |

## Columns (`*_hartree.csv`)

For the ortho–meta file, `isoA = ortho`, `isoB = meta`; for the para–meta
file, `isoA = para`, `isoB = meta`.

| Column | Meaning |
| --- | --- |
| `lambda_magnitude` | cavity coupling magnitude \|λ\| (a.u.) |
| `E_isoA_Hartrees`, `E_isoB_Hartrees` | converged electronic energy of isomer A / B (Hartree) |
| `zpe_isoA_Hartrees`, `zpe_isoB_Hartrees` | zero-point energy of A / B (Hartree) |
| `dE_isoA_isoB_raw_kcal/mol` | `E_A − E_B` (kcal/mol) |
| `dE_isoA_isoB_zpe_kcal/mol` | `(E_A + zpe_A) − (E_B + zpe_B)` (kcal/mol) |
| `isoA_converged`, `isoB_converged` | optimization convergence flag for each isomer |
| `isoA_has_freq`, `isoB_has_freq` | whether a frequency (ZPE) file exists for each isomer |

## Method notes

- Geometry type: **relaxed** — optimized under QED-DFT with the specific
  cavity direction and |λ| of that row (see `../../../README.md`).
- ZPE correction adds the harmonic zero-point energy to each isomer *before*
  differencing.
- All λ points are present and converged (see `availability.md`).
- Relative energies are pre-converted to kcal/mol (× 627.509) in the CSV;
  do **not** convert the `dE_*` columns again.

## Caveats

- The upstream `runs/` folder referenced by `analyze.py` is **not present**
  in `NEW_LAMBDA_SCAN/`; only these reduced outputs remain.
- `relative_energies_kcal.png` (written by `analyze.py`) is not currently
  generated; use `scripts/plot_isomer_energies.py` for figures.
