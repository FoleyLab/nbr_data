# outputs - Long QED-CCSD calculation outputs

This folder contains the original long QED-CCSD output files. They are
organized by geometry source and cavity direction:

| Folder | Pair | Direction | Notes |
| --- | --- | --- | --- |
| `relaxed_dir_70_31/` | ortho/meta | `(theta=70, phi=31)` | QED-DFT-relaxed geometries; lambda values 0.02-0.10. |
| `unrelaxed_dir_70_31/` | ortho/meta | `(theta=70, phi=31)` | Fixed unrelaxed geometries; lambda values 0.02-0.10. |
| `relaxed_dir_65_78/` | para/meta | `(theta=65, phi=78)` | QED-DFT-relaxed geometries; lambda values 0.02, 0.04, 0.08, 0.10. |
| `unrelaxed_dir_65_78/` | para/meta | `(theta=65, phi=78)` | Fixed unrelaxed geometries; lambda values 0.02-0.10. |

File names follow `<isomer>_lam<lambda>.out`. These files are calculation
provenance and should be treated as read-only source outputs. The reduced,
plot-ready data derived from them lives in `../summary/`.
