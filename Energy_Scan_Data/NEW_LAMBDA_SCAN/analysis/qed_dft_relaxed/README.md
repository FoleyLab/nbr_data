# qed_dft_relaxed — Relaxed QED-DFT energy tables

Per-lambda relative-energy tables for the two isomer pairs, computed from
geometries **relaxed under QED-DFT** for each (isomer, cavity direction,
|λ|) combination. Energies are in **Hartree**; only the per-point raw and
ZPE-corrected differences are given here (kcal/mol versions are produced by
`scripts/plot_isomer_energies.py`).

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

| Column | Meaning |
| --- | --- |
| `magnitude` | cavity coupling magnitude \|λ\| (a.u.) |
| `E_A_Ha`, `E_B_Ha` | converged electronic energy of isomer A / B (Hartree) |
| `zpe_A_Ha`, `zpe_B_Ha` | zero-point energy of A / B (Hartree) |
| `dE_raw_Ha` | `E_A − E_B` (Hartree) |
| `dE_zpe_Ha` | `(E_A + zpe_A) − (E_B + zpe_B)` (Hartree) |
| `A_converged`, `B_converged` | optimization convergence flag for each isomer |
| `A_has_freq`, `B_has_freq` | whether a frequency (ZPE) file exists for each isomer |

Isomer A is `ortho` / `para` for the two files respectively; B is `meta` in
both.

## Method notes

- Geometry type: **relaxed** — optimized under QED-DFT with the specific
  cavity direction and |λ| of that row (see `../../../README.md`).
- ZPE correction adds the harmonic zero-point energy to each isomer *before*
  differencing.
- All λ points are present and converged (see `availability.md`).
- To convert to kcal/mol: Hartree × 627.509.

## Caveats

- The upstream `runs/` folder referenced by `analyze.py` is **not present**
  in `NEW_LAMBDA_SCAN/`; only these reduced outputs remain.
- `relative_energies_kcal.png` (written by `analyze.py`) is not currently
  generated; use `scripts/plot_isomer_energies.py` for figures.
