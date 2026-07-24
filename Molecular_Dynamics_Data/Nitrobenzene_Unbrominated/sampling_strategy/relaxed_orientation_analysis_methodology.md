# Relaxed Orientation Energy Analysis Methodology

This note documents the sandbox analysis scripts used to compare sampled
QED-DFT relaxed intermediate energies against the older pQED unrelaxed
orientation surface.

## Input Data

The sampling reference is taken from this directory:

- `manifest.json`
- `frames.npz`

The sampled relaxed campaign summary is read from:

- `/Users/jfoley19/Code/nbr_data/Molecular_Dynamics_Data/Direction_A_MD_Sampled_Orientations/qed_dft_energy_summary.csv`

The complete sampled unrelaxed QED-DFT energies are read from:

- `/Users/jfoley19/Code/nbr_data/Molecular_Dynamics_Data/Direction_A_MD_Sampled_Orientations/UNRELAXED_Sampled_Orientations/direction_A_unrelaxed_results/direction_A_unrelaxed_energies.csv`

The pQED unrelaxed surface and reactant trajectory are read from:

- `/Users/jfoley19/Code/nbr_data/Molecular_Dynamics_Data/Nitrobenzene_Unbrominated/QED_DFT_wb97x_d/isomer_Nel_49_Nph_10_total_energies.dat`
- the `MD_FILE` currently configured in `/Users/jfoley19/Code/nbr_data/Molecular_Dynamics_Data/Nitrobenzene_Unbrominated/QED_DFT_wb97x_d/plot_timeseries_with_dwell_times.py`

At the time this note was written, `MD_FILE` points to the Direction A
trajectory, which is the correct trajectory for the sampled-orientation
campaign.

## Scripts

### `plot_orientation_energy_differences.py`

This script compares the sampled unrelaxed and relaxed QED-DFT energies on the
orientation support defined by `manifest.json` and `frames.npz`.

It:

1. Loads the 25 sampled orientations from `manifest.json`.
2. Loads the full folded trajectory support from `frames.npz`.
3. Reads the complete 75-cell unrelaxed QED-DFT CSV.
4. Reads the in-progress relaxed QED-DFT summary CSV.
5. Matches CSV rows to manifest orientations using the lambda/body-vector
   geometry rather than trusting rounded direction labels.
6. Computes electronic energy differences in kcal/mol:
   - `ortho - meta`
   - `para - meta`
7. For relaxed differences, uses only rows that are converged and have optimized
   XYZ files.
8. Writes orientation-map plots and a merged CSV.

The relaxed analysis in this script is conservative: stalled rows are excluded.

### `overlay_pqed_relaxed_timeseries.py`

This script overlays the sampled relaxed QED-DFT differences on the older pQED
unrelaxed orientation surface and on the trajectory time series used by
`plot_timeseries_with_dwell_times.py`.

It:

1. Reads the pQED grid from `isomer_Nel_49_Nph_10_total_energies.dat`.
2. Parses the trajectory selected by `MD_FILE` in
   `plot_timeseries_with_dwell_times.py`.
3. Applies the same theta branch convention used by that script:
   `theta > 100 -> 180 - theta`.
4. Reads all relaxed rows with `final_energy_hartree`, including stalled rows.
5. Computes provisional relaxed electronic energy differences:
   - `ortho - meta`
   - `para - meta`
6. Marks pairs as:
   - `converged` when both members are converged and optimized
   - `provisional_stalled` otherwise
7. Fits a low-order ridge-regularized Fourier/sinusoidal model to the available
   relaxed points.
8. Evaluates the fitted relaxed model along the selected trajectory.
9. Writes pQED surface overlays, fitted relaxed surfaces, a trajectory overlay,
   and CSV/JSON summaries.

This script is intentionally exploratory. Stalled rows are included so the
surface shape can be inspected before the campaign finishes. Those values should
be replaced by converged results as the campaign completes.

## Energy Convention

Both scripts use electronic energies only:

- relaxed: `final_energy_hartree`
- unrelaxed: `energy_hartree`
- pQED surface: intermediate total energies from the pQED grid

Energy differences are converted from Hartree to kcal/mol with:

```text
1 Hartree = 627.5094740631 kcal/mol
```

Zero-point energy is present in the relaxed summary file, but it is not included
in this pass. A later version can add parallel ZPE-corrected differences using
`zpe_corrected_energy_hartree`.

## Orientation Conventions

There are two coordinate conventions in play.

The sampling manifest uses body-frame orientation vectors `(a, b, c)` and folds
the site mirror with `b -> |b|`.

The older pQED time-series scripts use angular coordinates `(theta, phi)` and
apply a trajectory branch correction:

```text
theta > 100 degrees -> 180 - theta
```

For manifest-based orientation maps, the scripts use the folded body-frame
vectors from `frames.npz`.

For the pQED/time-series sandbox, sampled relaxed points are folded with the
same theta correction as `plot_timeseries_with_dwell_times.py` before overlaying
them on the pQED surface or trajectory.

## Outputs

`plot_orientation_energy_differences.py` writes:

- `unrelaxed_energy_differences_on_orientation_map.png`
- `relaxed_energy_differences_on_orientation_map.png`
- `relaxed_vs_unrelaxed_energy_differences.png`
- `orientation_energy_differences_merged.csv`
- `orientation_energy_differences_summary.json`

`overlay_pqed_relaxed_timeseries.py` writes:

- `pqed_surface_with_relaxed_samples.png`
- `relaxed_fourier_fit_surfaces.png`
- `trajectory_pqed_relaxed_timeseries_overlay.png`
- `relaxed_sandbox_pair_points.csv`
- `relaxed_points_nearest_trajectory_times.csv`
- `pqed_relaxed_overlay_summary.json`

## Current Caveats

- The relaxed campaign is incomplete.
- Meta is the limiting intermediate for paired relaxed differences.
- The Fourier/sinusoidal fit is useful as a smooth provisional guide, but it is
  not yet a final model.
- Stalled rows are useful in the sandbox overlay but should not be interpreted as
  final relaxed values.
- ZPE-corrected comparisons are a natural next pass once more frequencies are
  complete.
