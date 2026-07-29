# Audit Procedure for QED-CCSD Coordinates

## Audit Pass 1: Unrelaxed Geometries (COMPLETED)

**Status: All 20 files passed.** Every `meta`, `ortho`, and `para` `.out` file in both `unrelaxed_dir_70_31` and `unrelaxed_dir_65_78` uses coordinates identical to its source `.xyz` within the printed precision. No element or coordinate mismatches were found.

## Audit Pass 2: Relaxed Geometries (COMPLETED)

**Status: All 18 files passed.** Every `meta`, `para`, and `ortho` relaxed `.out` file uses coordinates identical to its corresponding source `.xyz`. No element or coordinate mismatches were found.

Perform the following audit to ensure that all relaxed QED-CCSD(2,2) calculations used consistent geometry coordinates from their respective source `.xyz` files.

### Audit Logic
For each `.out` file in a directory within `QED_CCSD/outputs/`, compare the lines under **"Geometry in angstrom"** with the corresponding `.xyz` file in the matching source path. The coordinates must be identical.

### Mapping convention
Source files are named `{isomer}_{dir}_{lam}.xyz` (e.g. `meta_65_78_lam0.02.xyz`). Output files are named `{isomer}_{lam}.out` (e.g. `meta_lam0.02.out`). The mapping is by isomer and lambda magnitude.

### Tasks

1.  **relaxed_dir_65_78 Meta** (4 files):
    -   `relaxed_dir_65_78/meta_65_78_lam0.{02,04,08,10}.xyz`
    -   → `QED_CCSD/outputs/relaxed_dir_65_78/meta_lam0.{02,04,08,10}.out`

2.  **relaxed_dir_65_78 Para** (4 files):
    -   `relaxed_dir_65_78/para_65_78_lam0.{02,04,08,10}.xyz`
    -   → `QED_CCSD/outputs/relaxed_dir_65_78/para_lam0.{02,04,08,10}.out`

3.  **relaxed_dir_70_31 Meta** (5 files):
    -   `relaxed_dir_70_31/meta_70_31_lam0.{02,04,06,08,10}.xyz`
    -   → `QED_CCSD/outputs/relaxed_dir_70_31/meta_lam0.{02,04,06,08,10}.out`

4.  **relaxed_dir_70_31 Ortho** (5 files):
    -   `relaxed_dir_70_31/ortho_70_31_lam0.{02,04,06,08,10}.xyz`
    -   → `QED_CCSD/outputs/relaxed_dir_70_31/ortho_lam0.{02,04,06,08,10}.out`

**Total: 18 files to verify.**

### Verification Step
For every file, confirm that the number of atoms and the spatial coordinates for each atom match perfectly between the output file's "Geometry in angstrom" block and the input `.xyz` file.