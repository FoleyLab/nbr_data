# Handoff Completeness and Accuracy Audit

Generated from `qed_ccsd_orientation_handoff.csv` and `trajectory_sampling_direction_A.dat`.

## Summary

- Expected entries: 75 total = 25 directions x 3 isomers (`ortho`, `meta`, `para`).
- CSV rows found: 50.
- Unique isomer/direction entries found: 50.
- Missing expected entries: 25.
- Redundant isomer/direction entries: 0.
- Unexpected isomer/direction entries: 0.
- Rows with any accuracy issue: 0.
- Referenced `.xyz` files missing from Directory 1: 0.
- Source `optimized.xyz` files missing from `runs_trajectory`: 0.
- Coordinate comparisons completed: 50.
- Coordinate mismatches vs source `optimized.xyz`: 0.
- Full-file mismatches vs source `optimized.xyz`: 0.

## Completeness

| Isomer | Present count | Missing count | Present directions | Missing directions |
|---|---:|---:|---|---|
| ortho | 25 | 0 | A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19, A20, A21, A22, A23, A24, A25 | None |
| meta | 4 | 21 | A2, A6, A8, A11 | A1, A3, A4, A5, A7, A9, A10, A12, A13, A14, A15, A16, A17, A18, A19, A20, A21, A22, A23, A24, A25 |
| para | 21 | 4 | A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18, A21, A22, A25 | A19, A20, A23, A24 |

### Missing Entries

| Isomer | Direction | Expected theta | Expected phi | Expected lambda_magnitude | Expected xyz_file_name |
|---|---|---:|---:|---:|---|
| meta | A1 | 101.0 | 31.0 | 0.1 | `meta_dirA1_r001_th101_ph31.xyz` |
| meta | A3 | 113.0 | 36.0 | 0.1 | `meta_dirA3_r001_th113_ph36.xyz` |
| meta | A4 | 125.0 | 29.0 | 0.1 | `meta_dirA4_r001_th125_ph29.xyz` |
| meta | A5 | 136.0 | 36.0 | 0.1 | `meta_dirA5_r001_th136_ph36.xyz` |
| meta | A7 | 125.0 | 50.0 | 0.1 | `meta_dirA7_r001_th125_ph50.xyz` |
| meta | A9 | 138.0 | 55.0 | 0.1 | `meta_dirA9_r001_th138_ph55.xyz` |
| meta | A10 | 124.0 | 69.0 | 0.1 | `meta_dirA10_r001_th124_ph69.xyz` |
| meta | A12 | 149.0 | 65.0 | 0.1 | `meta_dirA12_r001_th149_ph65.xyz` |
| meta | A13 | 139.0 | 83.0 | 0.1 | `meta_dirA13_r001_th139_ph83.xyz` |
| meta | A14 | 157.0 | 73.0 | 0.1 | `meta_dirA14_r001_th157_ph73.xyz` |
| meta | A15 | 134.0 | 71.0 | 0.1 | `meta_dirA15_r001_th134_ph71.xyz` |
| meta | A16 | 159.0 | 40.0 | 0.1 | `meta_dirA16_r001_th159_ph40.xyz` |
| meta | A17 | 146.0 | 42.0 | 0.1 | `meta_dirA17_r001_th146_ph42.xyz` |
| meta | A18 | 153.0 | 13.0 | 0.1 | `meta_dirA18_r001_th153_ph13.xyz` |
| meta | A19 | 130.0 | 39.0 | 0.1 | `meta_dirA19_r001_th130_ph39.xyz` |
| meta | A20 | 139.0 | 1.0 | 0.1 | `meta_dirA20_r001_th139_ph1.xyz` |
| meta | A21 | 135.0 | 21.0 | 0.1 | `meta_dirA21_r001_th135_ph21.xyz` |
| meta | A22 | 116.0 | 41.0 | 0.1 | `meta_dirA22_r001_th116_ph41.xyz` |
| meta | A23 | 119.0 | 26.0 | 0.1 | `meta_dirA23_r001_th119_ph26.xyz` |
| meta | A24 | 124.0 | 9.0 | 0.1 | `meta_dirA24_r001_th124_ph9.xyz` |
| meta | A25 | 105.0 | 32.0 | 0.1 | `meta_dirA25_r001_th105_ph32.xyz` |
| para | A19 | 130.0 | 39.0 | 0.1 | `para_dirA19_r001_th130_ph39.xyz` |
| para | A20 | 139.0 | 1.0 | 0.1 | `para_dirA20_r001_th139_ph1.xyz` |
| para | A23 | 119.0 | 26.0 | 0.1 | `para_dirA23_r001_th119_ph26.xyz` |
| para | A24 | 124.0 | 9.0 | 0.1 | `para_dirA24_r001_th124_ph9.xyz` |

No redundant isomer/direction entries were found.

## Accuracy Checks

- Spherical-polar to Cartesian check tolerance: `1e-12`.
- Maximum component residual across all rows: `1.388e-17`.
- Maximum norm residual, `sqrt(Ex**2 + Ey**2 + Ez**2) - lambda_magnitude`: `2.776e-17`.
- Maximum residual from required lambda magnitude `0.1`: `0.000e+00`.
- Filename theta/phi checks passed for rows without listed issues below.
- Coordinate block comparisons passed for 50 of 50 checked rows.

No row-level accuracy issues were found.

## Interpretation

The handoff CSV is not yet complete under the requested 75-entry criterion. The accuracy checks on the rows that are present are otherwise summarized above; any row-level problems are listed explicitly.

For the remaining related goals, a new chat would make sense if they involve a different output artifact or a broader analysis/plotting workflow. If they build directly on this audit file and the same CSVs, continuing here is also fine.
