# Projected Gradient Ring Analysis

`analyze_projected_gradient_ring.py` tests whether localized jagged features in
orientation-dependent energy-difference surfaces correlate with orientations
where the projected-out torque/constraint component becomes small for one
intermediate.

The script reads optimized constrained geometries from folders like:

```text
ortho_th72_ph180/optimized.xyz
ortho_th72_ph180/cell.json
meta_th72_ph180/optimized.xyz
meta_th72_ph180/cell.json
```

For each requested `(theta, phi, isomer)` case, it:

1. Parses `optimized.xyz` into a Psi4 geometry string.
2. Reads `cell.json["lambda_vector"]`.
3. Runs CQED-DFT full-gradient and projected-gradient calculations.
4. Computes:

```python
norm_removed = np.linalg.norm(g_full - g_projected)
```

This `norm_removed` value measures the gradient component removed by the
projection. If a species naturally has vanishing torque at a given orientation,
`norm_removed` should become locally small for that species.

## Dry Run

Use a dry run first to verify folder discovery and lambda-vector parsing without
running Psi4:

```bash
cd /Users/jfoley19/Code/nbr_data/QED_DFT_Theta_Phi_Relax_Campaign/runs_grid_no_freq

/Users/jfoley19/miniforge3/bin/conda run -n p4dev python analyze_projected_gradient_ring.py \
  --thetas 72 \
  --phis 180:270:18 \
  --isomers ortho meta \
  --dry-run
```

## Run The Analysis

```bash
cd /Users/jfoley19/Code/nbr_data/QED_DFT_Theta_Phi_Relax_Campaign/runs_grid_no_freq

MPLCONFIGDIR=/private/tmp /Users/jfoley19/miniforge3/bin/conda run -n p4dev python analyze_projected_gradient_ring.py \
  --thetas 72 \
  --phis 180:270:18 \
  --isomers ortho meta \
  --plot
```

The `MPLCONFIGDIR=/private/tmp` setting avoids Matplotlib cache-permission
warnings when writing the optional diagnostic plot.

## Useful Variants

Scan multiple theta rings:

```bash
MPLCONFIGDIR=/private/tmp /Users/jfoley19/miniforge3/bin/conda run -n p4dev python analyze_projected_gradient_ring.py \
  --thetas 63:81:9 \
  --phis 180:270:18 \
  --isomers ortho meta \
  --output-prefix projected_gradient_theta63_81_phi180_270 \
  --plot
```

Use explicit phi values:

```bash
MPLCONFIGDIR=/private/tmp /Users/jfoley19/miniforge3/bin/conda run -n p4dev python analyze_projected_gradient_ring.py \
  --thetas 72 \
  --phis 180,198,216,234,252,270 \
  --isomers ortho meta \
  --plot
```

Compare a different pair:

```bash
MPLCONFIGDIR=/private/tmp /Users/jfoley19/miniforge3/bin/conda run -n p4dev python analyze_projected_gradient_ring.py \
  --thetas 72 \
  --phis 180:270:18 \
  --isomers para meta \
  --output-prefix projected_gradient_para_meta_theta72 \
  --plot
```

## Outputs

For the default `--output-prefix projected_gradient_ring`, the script writes:

```text
projected_gradient_ring_raw.csv
projected_gradient_ring_ortho_meta_pairwise.csv
projected_gradient_ring_ortho_meta_diagnostic.png
```

The raw CSV has one row per `(theta, phi, isomer)`:

```text
theta
phi
isomer
energy_full_hartree
energy_projected_hartree
delta_energy_projected_minus_full_hartree
norm_full
norm_projected
norm_removed
lambda_x
lambda_y
lambda_z
geometry_path
```

The pairwise CSV compares the first two isomers passed to `--isomers`, including:

```text
delta_energy_a_minus_b_kcal_mol
norm_removed_a
norm_removed_b
delta_norm_removed_a_minus_b
```

For the ortho-meta hypothesis:

- A blue feature in `E_ortho - E_meta` should coincide with a local dip in
  `norm_removed_ortho` if the ortho torque vanishes there.
- A red feature in `E_ortho - E_meta` should coincide with a local dip in
  `norm_removed_meta` if the meta torque vanishes there.

## Relevant Options

```text
--root              runs_grid_no_freq folder; defaults to this script's folder
--thetas            comma list or start:stop:step, e.g. 72 or 63:81:9
--phis              comma list or start:stop:step, e.g. 180:270:18
--isomers           isomer names; first two are used for pairwise output
--output-prefix     prefix for output CSV/PNG files
--dry-run           validate inputs without Psi4 calculations
--plot              write the quick pairwise diagnostic PNG
--memory            Psi4 memory, default: 4 GB
--num-threads       Psi4 thread count, default: 1
```
