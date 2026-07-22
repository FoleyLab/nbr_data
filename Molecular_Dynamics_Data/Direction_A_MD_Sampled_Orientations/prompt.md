# Recreate QED-DFT Orientation Data Handoff

You can work in the following directories:

Directory 1, read access: /Users/jfoley19/Code/qcage/nitrobenzene_opt_and_freq/runs_trajectory

Directory 2, write access: /Users/jfoley19/Code/nbr_data/Molecular_Dynamics_Data/Direction_A_MD_Sampled_Orientations

Directory 1 contains subdirectories with QED-DFT geometry optimization and frequency results for para and meta bromonitrobenzene orientations. The subdirectories are named like:

`para_dirA13_r001_th139_ph83` for different directions for para, `meta_dirA10_r001_th124_ph69` for different directions for meta.

Each subdirectory should contain `cell.json` and `opt_status.json`; some contain `optimized.xyz`; some contain frequency products such as `frequencies.json`, `freq_cm.npy`, `H_cart_au.npy`, and `H_locked_mw.npy`.

Please compile updated data products in Directory 2.

## Direction Validation

Before writing the CSV files, verify that each `cell.json` maps `theta`, `phi`, and `magnitude` to `lambda_vector` using the spherical-polar to Cartesian transformation, with angles in degrees:

`Ex = magnitude * sin(theta) * cos(phi)`

`Ey = magnitude * sin(theta) * sin(phi)`

`Ez = magnitude * cos(theta)`

Also verify:

`sqrt(Ex**2 + Ey**2 + Ez**2) = magnitude`

Treat differences near floating-point precision as acceptable, but report any run with a meaningful mismatch.

## Goal 1: QED-CCSD Geometry Handoff

Create `qed_ccsd_orientation_handoff.csv` with one row for each subdirectory that contains `optimized.xyz`.

Columns:

`theta,phi,Ex,Ey,Ez,lambda_magnitude,xyz_file_name`

For each row:

- `isomer` comes from `cell.json["isomer"]`.
- `theta` comes from `cell.json["theta"]`.
- `phi` comes from `cell.json["phi"]`.
- `Ex`, `Ey`, and `Ez` come from `cell.json["lambda_vector"][0]`, `[1]`, and `[2]`.
- `lambda_magnitude` comes from `cell.json["magnitude"]`.
- `xyz_file_name` should be `{cell.json["id"]}.xyz`.

Also copy each source `optimized.xyz` into Directory 2 and rename the copy to `{cell.json["id"]}.xyz`. Preserve the contents of `optimized.xyz` verbatim, including the second-line comment.

## Goal 2: QED-DFT Plotting Summary

Create `qed_dft_energy_summary.csv` with one row for every subdirectory in Directory 1, including runs without `optimized.xyz`.

Include columns useful for plotting energy vs configuration with and without zero point correction:

`isomer,id,direction_label,theta,phi,Ex,Ey,Ez,lambda_magnitude,has_optimized_xyz,converged,promoted,resumed_from_stall,attempts,final_gnorm,conv_threshold,final_energy_hartree,frequency_complete,zpe_hartree,zpe_ev,zpe_corrected_energy_hartree,n_real_modes,n_imaginary,imaginary_freqs_cm,lowest_freq_cm,xyz_file_name`

Use `opt_status.json` for optimization status and final energy. Use `frequencies.json` for ZPE and imaginary-frequency fields when present. If no frequency data exists, leave ZPE/frequency fields blank. Compute:

`zpe_corrected_energy_hartree = final_energy_hartree + zpe_hartree`

when both values are available.

Before finishing, verify:

- The spherical-polar to Cartesian direction check passes for every `cell.json`.
- The norm of each `lambda_vector` matches `cell.json["magnitude"]`.
- The handoff CSV row count equals the number of copied `.xyz` files.
- The energy summary row count equals the number of top-level run subdirectories.
- Every `xyz_file_name` in the handoff CSV exists in Directory 2.
- No source files in Directory 1 were modified.
