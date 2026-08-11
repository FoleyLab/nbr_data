# NEW_LAMBDA_SCAN Dataset

This directory contains geometry inputs and organized energy results for
cavity-coupled calculations on three bromonitrobenzene intermediates:
ortho-, meta-, and para-bromonitrobenzene.

The study focuses on relative isomer energetics under strong directional
cavity coupling. Two cavity directions are used:

| Direction | Primary comparison | Geometry/energy pair |
| --- | --- | --- |
| `(theta=70, phi=31)` | `ortho - meta` | ortho and meta |
| `(theta=65, phi=78)` | `para - meta` | para and meta |

## Common Calculation Configuration

| Setting | Value |
| --- | --- |
| Basis set | `6-311G*` |
| Cavity frequency | `0.06615` a.u. |
| Charge | `+1` |
| Multiplicity | `1` |
| Main correlated method | QED-CCSD(2,2) |
| Geometry relaxation method | QED-DFT |

## Directory Layout

```text
NEW_LAMBDA_SCAN/
├── README.md
├── GEOMETRIES/
│   ├── relaxed_dir_70_31/
│   ├── relaxed_dir_65_78/
│   ├── unrelaxed_dir_70_31/
│   └── unrelaxed_dir_65_78/
└── Lambda_Scan_Results/
    ├── README.md
    ├── qed_dft_relaxed/
    ├── qed_dft_unrelaxed/
    ├── QED_CCSD/
    ├── pqed/
    ├── scripts/
    └── docs/
```

## `GEOMETRIES/`

`GEOMETRIES/` contains the `.xyz` structures and `lambda_grid.csv` files used
as inputs for relaxed and unrelaxed calculations.

| Folder | Geometry type | Pair | Direction | Lambda points |
| --- | --- | --- | --- | --- |
| `GEOMETRIES/unrelaxed_dir_70_31/` | unrelaxed | ortho/meta | `(70,31)` | 0.02, 0.04, 0.06, 0.08, 0.10 |
| `GEOMETRIES/unrelaxed_dir_65_78/` | unrelaxed | para/meta | `(65,78)` | 0.02, 0.04, 0.06, 0.08, 0.10 |
| `GEOMETRIES/relaxed_dir_70_31/` | QED-DFT relaxed | ortho/meta | `(70,31)` | 0.02, 0.04, 0.06, 0.08, 0.10 |
| `GEOMETRIES/relaxed_dir_65_78/` | QED-DFT relaxed | para/meta | `(65,78)` | 0.02, 0.04, 0.08, 0.10 |

### Unrelaxed Geometries

The unrelaxed folders contain fixed reference geometries:

```text
GEOMETRIES/unrelaxed_dir_70_31/
  lambda_grid.csv
  ortho.xyz
  meta.xyz

GEOMETRIES/unrelaxed_dir_65_78/
  lambda_grid.csv
  para.xyz
  meta.xyz
```

For unrelaxed scans, each `.xyz` file is reused for every row of the matching
`lambda_grid.csv`.

### Relaxed Geometries

The relaxed folders contain QED-DFT-optimized geometries for specific
`(isomer, direction, |lambda|)` points:

```text
GEOMETRIES/relaxed_dir_70_31/
  lambda_grid.csv
  ortho_70_31_lam0.02.xyz ... ortho_70_31_lam0.10.xyz
  meta_70_31_lam0.02.xyz  ... meta_70_31_lam0.10.xyz

GEOMETRIES/relaxed_dir_65_78/
  lambda_grid.csv
  para_65_78_lam0.02.xyz ... para_65_78_lam0.10.xyz
  meta_65_78_lam0.02.xyz ... meta_65_78_lam0.10.xyz
```

`relaxed_dir_65_78` does not contain `lam0.06` relaxed geometries; that absence
is intentional in this dataset.

The second line of each relaxed `.xyz` records the isomer, direction, lambda
magnitude, optimization status, gradient information, and QED-DFT energy.
Some geometries are marked `PROMOTED` rather than `OPTIMIZED`; those were
accepted because their gradient norms are below `5.0e-04`.

## Lambda Grid Files

Each geometry subfolder has a `lambda_grid.csv` with columns:

```text
theta,phi,Ex,Ey,Ez,lambda_magnitude,file_A,file_B
```

| Column | Meaning |
| --- | --- |
| `theta`, `phi` | cavity coupling-vector direction in degrees |
| `Ex`, `Ey`, `Ez` | Cartesian components of the lambda vector |
| `lambda_magnitude` | vector magnitude `|lambda|` in a.u. |
| `file_A`, `file_B` | geometry files to use for that row |

The Cartesian vector components are generated from:

```python
Ex = lambda_magnitude * sin(theta) * cos(phi)
Ey = lambda_magnitude * sin(theta) * sin(phi)
Ez = lambda_magnitude * cos(theta)
```

with `theta` and `phi` interpreted in degrees.

## `Lambda_Scan_Results/`

`Lambda_Scan_Results/` contains the organized energy outputs and plotting
tools derived from these geometries.

| Folder | Contents |
| --- | --- |
| `qed_dft_relaxed/` | QED-DFT energies at the QED-DFT-relaxed geometries, including ZPE where available. |
| `qed_dft_unrelaxed/` | QED-DFT energies at fixed unrelaxed geometries. |
| `QED_CCSD/` | QED-CCSD(2,2), QED-HF, and long original calculation outputs. |
| `pqed/` | pQED unrelaxed scans for 49 electrons, `Nph=3/10`, with and without coherent-state transform. |
| `scripts/` | Plotting script for method and geometry comparisons. |
| `docs/` | Presentation material derived from the summarized results. |

See `Lambda_Scan_Results/README.md` for the detailed energy-data
organization, CSV column conventions, and plotting workflow.

## Relationship Between Geometries and Energies

- `GEOMETRIES/relaxed_*` provides the QED-DFT-relaxed structures used for
  relaxed QED-DFT summaries and relaxed-geometry QED-CCSD calculations.
- `GEOMETRIES/unrelaxed_*` provides the fixed reference structures used for
  unrelaxed QED-DFT, unrelaxed QED-CCSD, and pQED scans.
- pQED is available only at unrelaxed geometries.
- Long QED-CCSD calculation outputs are preserved in
  `Lambda_Scan_Results/QED_CCSD/outputs/`; reduced plot-ready summaries are in
  `Lambda_Scan_Results/QED_CCSD/summary/`.

## Pre-Handoff Consistency Checks

The geometry inputs were checked for the following before handoff:

- Each geometry folder contains the expected isomer subset.
- Each `lambda_grid.csv` direction matches the `dir_<theta>_<phi>` string in
  the folder name.
- The `Ex`, `Ey`, and `Ez` values are consistent with the listed `theta`,
  `phi`, and `lambda_magnitude`.
- The norm of each lambda vector matches `lambda_magnitude`.
- Every `file_A` and `file_B` entry names an existing `.xyz` file in the same
  geometry folder.
- Relaxed `.xyz` filenames match the direction encoded in their folder.
- Relaxed `.xyz` comment lines match the corresponding filename stem.
- All `.xyz` files contain 15 atoms with composition `C6 H5 Br N O2`.
