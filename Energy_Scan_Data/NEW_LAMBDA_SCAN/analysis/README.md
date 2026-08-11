# analysis — Reorganized results for the cavity opt+frequency campaign

Relative energies of bromonitrobenzene isomer pairs vs. cavity-coupling
magnitude |λ|, along one cavity orientation per pair, from multiple levels of
theory. Reorganized by method (see `git log`/subfolders); each subfolder has
its own README with a file-level table.

## Layout

```text
analysis/
├── README.md                 # this file (index + method overview)
├── qed_dft_relaxed/          # relaxed-geometry QED-DFT tables (Hartree), from analyze.py
├── qed_dft_unrelaxed/        # unrelaxed-geometry QED-DFT scans (kcal/mol)
├── pqed/                     # parameterized QED scans, Nel=49, Nph=3/10, with/without CS
├── polariton_scan/           # full-sphere polariton orientation scans (DIFFERENT SYSTEM)
├── scripts/                  # plotting/analysis tools
└── docs/                     # presentations
```

## Files by folder (quick index)

| Folder | Files | Units | Produced by |
| --- | --- | --- | --- |
| `qed_dft_relaxed/` | `ortho-meta_70_31_hartree.csv`, `para-meta_65_78_hartree.csv`, `availability.md` | Hartree | `../../analyze.py` (from `runs/`) |
| `qed_dft_unrelaxed/` | `ortho_meta_dir70_31_scan_qed_dft_no_relax.csv`, `para_meta_dir65_78_scan_qed_dft_no_relax.csv` | Hartree / kcal/mol | external QED-DFT scans |
| `pqed/` | 8 × `pqed_49_{3,10}_dir{70_31,65_78}_scan[_CS].csv` | Hartree / kcal/mol | pQED scans (also in `../../QED_CCSD/summary/`) |
| `polariton_scan/` | `polariton_energy_scan.csv`, `polariton_energy_scan_fc.csv` | Hartree | unknown (see README) |
| `scripts/` | `plot_isomer_energies.py` + its 4 generated figures | — | — |
| `docs/` | `relaxed_vs_unrelaxed.pptx` | — | — |

## Method overview

- **System:** three bromonitrobenzene isomers (ortho-, meta-, para-), C6H4BrNO2.
  Basis 6-311G\*, charge +1, singlet. Level of theory for the relaxed/unrelaxed
  CC energies: QED-CCSD(2,2); relaxation geometries from QED-DFT.
- **Cavity directions:** (θ=70°, φ=31°) for the ortho–meta pair and
  (θ=65°, φ=78°) for the para–meta pair. λ vector from spherical coordinates:
  `Ex = |λ|·sinθ·cosφ`, `Ey = |λ|·sinθ·sinφ`, `Ez = |λ|·cosθ`.
- **λ magnitudes:** 0.02, 0.04, 0.06, 0.08, 0.10 a.u. (relaxed dir_65_78 is
  missing λ=0.06 by design — see `../../README.md`).
- **Relative energies:** ΔE = E(A) − E(B), reported for `ortho − meta` and
  `para − meta`. `qed_dft_relaxed/` additionally has ZPE-corrected differences,
  `(E_A + zpe_A) − (E_B + zpe_B)`.
- **Units:** Hartree tables live in `qed_dft_relaxed/`; everything in
  kcal/mol uses 627.509 Ha→kcal/mol.
- **Relaxed vs unrelaxed:** relaxed geometries are re-optimized under QED-DFT
  per (isomer, direction, |λ|); unrelaxed use a fixed gas-phase geometry at
  every |λ|.

## Relationships to the parent folder

| Parent item | Role |
| --- | --- |
| `analyze.py` | regenerates `qed_dft_relaxed/*` from `runs/` (folder currently absent) |
| `plot_energies.py` | QED-CCSD (relaxed/unrelaxed) vs pQED figures from `QED_CCSD/summary/` |
| `plot_qed_ccsd_energies_relaxed_vs_unrelaxed.py` | relaxed-vs-unrelaxed QED-CCSD figures |
| `QED_CCSD/summary/` | canonical pQED CSVs + QED-CCSD/QEDHF energy tables |
| `README.md`, `audit_procedure.md` | dataset documentation and geometry audits |

## Known gaps / caveats

1. `qed_dft_relaxed/` data came from a `runs/` folder that is no longer present;
   only the reduced tables survive.
2. `pqed/` duplicates `QED_CCSD/summary/unrelaxed_dir_*_pqed_49_*.csv` — keep in
   sync (see `pqed/README.md`).
3. `polariton_scan/` belongs to a different molecule (~ −435 Ha) and has no
   generating script here — see its README before use.
