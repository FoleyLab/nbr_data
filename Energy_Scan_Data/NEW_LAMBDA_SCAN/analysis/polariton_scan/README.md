# polariton_scan — Full spherical polariton energy scans

Two near-identical full-sphere scans of the polariton energy vs. cavity
orientation. Each has **8100 rows = 90 θ × 90 φ** (θ 0–180°, φ 0–360°),
i.e. one point per ~2° step in each angle.

## Files

| File | Notes |
| --- | --- |
| `polariton_energy_scan.csv` | base variant |
| `polariton_energy_scan_fc.csv` | `fc` variant; differs from the base file in **every** value |

## Columns

| Column | Meaning |
| --- | --- |
| `theta`, `phi` | cavity direction (degrees) |
| `E_el` | electronic energy in cavity (depends on orientation) |
| `E_el_uncoupled` | gas-phase (uncoupled) electronic reference; a single constant across all rows |
| `E_ph` | photon energy contribution |
| `E_blc` | bilinear coupling contribution |
| `E_dse` | dipole self-energy contribution |
| `E_total` | total polariton energy (≈ `E_el + E_ph + E_blc + E_dse`) |

## Caveats (important)

- **Different molecule/system.** Energies are ~ −435 Ha, whereas the
  bromonitrobenzene data elsewhere in this `analysis/` folder is ~ −3010 Ha.
  These files do not belong to the nitrobenzene isomer analysis.
- **Unknown provenance.** No generating script for these files exists in this
  directory or in `NEW_LAMBDA_SCAN/`. The meaning of the `fc` suffix has not
  been verified.
- Handle/interpret with care; do not mix with the isomer-relative-energy data.
