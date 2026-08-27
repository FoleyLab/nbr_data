# NEW_LAMBDA_SCAN Dataset

This dataset contains geometry inputs and organized energy results for cavity-coupled calculations on three bromonitrobenzene intermediates (ortho-, meta-, and para-). The study identifies the effects of strong directional coupling on relative isomer energetics.

## Quick Start
To generate publication-quality plots from the prepared summary data:
```bash
python plot_qed_ccsd_relaxed_unrelaxed_styled.py
```
This will produce `qed_ccsd_relaxed_vs_unrelaxed_styled.png` and `.pdf` files in the root directory.

## Directory Overview
- **GEOMETRIES/**: Input structures and grid definitions.
  - `relaxed_dir_...`: QED-DFT optimized geometries for specific coupling directions.
  - `unrelaxed_dir_...`: Fixed reference gas-phase geometries used for baseline comparisons.
- **Lambda_Scan_Results/**: Processed results and analysis tools.
  - `QED_CCSD/`: Primary results including QED-CCSD(2,2) and QED-HF calculations.
  - `qed_dft_unrelaxed/`: Results for the baseline unrelaxed scans.
  - `scripts/`: Python scripts for data visualization and processing.
  - `docs/`: Documentation materials for publication prep.

## Methods
Calculations were performed using the following configurations:
- **Basis Set**: 6-311G*
- **Cavity Frequency**: 0.06615 a.u.
- **Charge/Multiplicity**: +1 / 1
- **Primary Correlated Method**: QED-CCSD(2,2)
- **Geometry Optimization**: QED-DFT for relaxed structures.

## File Formats
- **.xyz**: Standard coordinate files. For "relaxed" samples, the second line contains metadata (isomer, direction, $\lambda$, status).
- **.csv**: Machine-readable tables. 
  - Geometry folders contain `lambda_grid.csv` defining $(\theta, \phi)$ and magnitude of the coupling vector.
  - Result directories contain summary files with standardized columns: magnitude, results in Hartrees/kcal/mol, and calculated $\Delta E$.

## Reproducibility
Data integrity is ensured by several consistency checks:
- **Grid Consistency**: Every entry's $(\theta, \phi)$ coordinates exactly match the folder naming convention.
- **Vector Normalization**: The Cartesian components $(E_x, E_y, E_z)$ are consistent with the `lambda_magnitude`.
- **Conversion Constants**: 1 unit of $\Delta E$ in the raw data is consistent between Hartrees and kcal/mol conversions (627.509).

## Validation
Pre-handoff checks verify:
1. Every `lambda_grid.csv` contains all expected isomer entries.
2. Reference geometry files exist for every point on the grid.
3. All `.xyz` files contain exactly 15 atoms ($C_6 H_5 Br N O_2$).
4. The summary results match the dimensions of the input `lambda_grid` files.
