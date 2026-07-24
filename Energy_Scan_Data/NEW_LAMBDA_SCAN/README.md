# NEW_LAMBDA_SCAN Dataset

This directory contains geometry and optical-cavity coupling-vector inputs for
CQED coupled-cluster calculations on three bromonitrobenzene isomers:
ortho-, meta-, and para-bromonitrobenzene.

The dataset is organized by geometry type and by the direction of the cavity
coupling vector. Each row of `lambda_grid.csv` gives the cavity vector and the
two geometry files to use for that calculation.

## Common configurations
Basis set: 6-311G*
frequency: 0.06615 atomic units
Level of theory: QED-CCSD(2,2)
Charge: +1
Multiplicity: 1

## Directory Contents

```text
NEW_LAMBDA_SCAN/
  unrelaxed_dir_70_31/
    lambda_grid.csv
    ortho.xyz
    meta.xyz

  unrelaxed_dir_65_78/
    lambda_grid.csv
    meta.xyz
    para.xyz

  relaxed_dir_70_31/
    lambda_grid.csv
    ortho_70_31_lam0.02.xyz
    ortho_70_31_lam0.04.xyz
    ortho_70_31_lam0.06.xyz
    ortho_70_31_lam0.08.xyz
    ortho_70_31_lam0.10.xyz
    meta_70_31_lam0.02.xyz
    meta_70_31_lam0.04.xyz
    meta_70_31_lam0.06.xyz
    meta_70_31_lam0.08.xyz
    meta_70_31_lam0.10.xyz

  relaxed_dir_65_78/
    lambda_grid.csv
    meta_65_78_lam0.02.xyz
    meta_65_78_lam0.04.xyz
    meta_65_78_lam0.08.xyz
    meta_65_78_lam0.10.xyz
    para_65_78_lam0.02.xyz
    para_65_78_lam0.04.xyz
    para_65_78_lam0.08.xyz
    para_65_78_lam0.10.xyz
```

## Geometry Types

### Unrelaxed Geometries

The `unrelaxed_*` directories contain fixed gas-phase/unrelaxed geometries.
For these cases, use each `.xyz` geometry with every row of the corresponding
`lambda_grid.csv` file in the same directory.

- `unrelaxed_dir_70_31` contains ortho and meta geometries.
- `unrelaxed_dir_65_78` contains meta and para geometries.

### Relaxed Geometries

The `relaxed_*` directories contain geometries optimized under CQED-DFT with a
specific cavity coupling-vector direction and magnitude. In these directories,
each `.xyz` file already corresponds to one specific isomer, direction, and
lambda magnitude.

For relaxed calculations, use each relaxed `.xyz` file with the matching row in
the `lambda_grid.csv` file from the same directory. The lambda magnitude is
encoded in the filename as `lam0.XX`.

- `relaxed_dir_70_31` contains ortho and meta relaxed geometries.
- `relaxed_dir_65_78` contains meta and para relaxed geometries.
- The absence of `lam0.06` files in `relaxed_dir_65_78` is intentional.

## Lambda Grid Files

Each `lambda_grid.csv` has the following columns:

```text
theta,phi,Ex,Ey,Ez,lambda_magnitude,file_A,file_B
```

where:

- `theta` and `phi` define the cavity coupling-vector direction in degrees.
- `Ex`, `Ey`, and `Ez` are the Cartesian components of the lambda vector.
- `lambda_magnitude` is the vector magnitude.
- `file_A` and `file_B` name the two `.xyz` geometry files to use for that row.

The directory name encodes the intended direction:

- `*_dir_70_31` means `theta = 70`, `phi = 31`.
- `*_dir_65_78` means `theta = 65`, `phi = 78`.

The Cartesian components are generated from the spherical direction and scaled
to the requested magnitude:

```python
Ex = lambda_magnitude * sin(theta) * cos(phi)
Ey = lambda_magnitude * sin(theta) * sin(phi)
Ez = lambda_magnitude * cos(theta)
```

with `theta` and `phi` interpreted in degrees.

All lambda magnitudes are in the range `0.02` to `0.10`.

## Intended Workflow

### 1. Unrelaxed Scans

For each `unrelaxed_*` directory:

1. Read the direction and lambda vectors from `lambda_grid.csv`.
2. For each row, use the unrelaxed geometries named in `file_A` and `file_B`.
3. Run the intended CQED coupled-cluster calculation for the listed lambda
   vector and geometry pair.
4. Keep the output labeled by isomer pair, direction, and lambda magnitude.

Expected unrelaxed combinations:

```text
unrelaxed_dir_70_31:
  file_A = ortho.xyz
  file_B = meta.xyz
  lambda = 0.02, 0.04, 0.06, 0.08, 0.10

unrelaxed_dir_65_78:
  file_A = para.xyz
  file_B = meta.xyz
  lambda = 0.02, 0.04, 0.06, 0.08, 0.10
```

### 2. Relaxed-Geometry Calculations

For each `relaxed_*` directory:

1. Read the direction and lambda vectors from `lambda_grid.csv`.
2. For each row, use the relaxed geometries named in `file_A` and `file_B`.
3. Run the intended CQED coupled-cluster calculation for the listed lambda
   vector and geometry pair.

Expected relaxed combinations:

```text
relaxed_dir_70_31:
  file_A = ortho_70_31_lam0.XX.xyz
  file_B = meta_70_31_lam0.XX.xyz
  lambda = 0.02, 0.04, 0.06, 0.08, 0.10

relaxed_dir_65_78:
  file_A = para_65_78_lam0.XX.xyz
  file_B = meta_65_78_lam0.XX.xyz
  lambda = 0.02, 0.04, 0.08, 0.10
```

## Notes on Relaxed XYZ Comments

The second line of each relaxed `.xyz` file records the isomer, direction,
lambda magnitude, optimization status, gradient information, and energy.

Some geometries are marked `PROMOTED` rather than `OPTIMIZED`. These were
accepted because their gradient norms are below `5.0e-04`, which is sufficient
for the intended downstream calculations.

## Pre-Handoff Consistency Checks

The dataset was checked for the following before handoff:

- Each directory contains the expected isomer subset.
- Each `lambda_grid.csv` direction matches the `dir_{theta}_{phi}` string in
  the directory name.
- The `Ex`, `Ey`, and `Ez` values are consistent with the listed `theta`,
  `phi`, and `lambda_magnitude`.
- The norm of each lambda vector matches `lambda_magnitude`.
- Every `file_A` and `file_B` entry names an existing `.xyz` file in the same
  directory.
- Relaxed `.xyz` filenames match the direction encoded in their directory.
- Relaxed `.xyz` comment lines match the corresponding filename stem.
- All `.xyz` files contain 15 atoms with composition `C6 H5 Br N O2`.
- No editor swap files or common temporary files are present.

