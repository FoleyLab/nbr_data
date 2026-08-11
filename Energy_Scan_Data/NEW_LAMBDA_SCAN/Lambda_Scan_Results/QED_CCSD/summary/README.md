# summary - Reduced QED-CCSD, QED-HF, and pQED tables

This folder contains compact summary tables used for plotting and method
comparisons.

## QED-CCSD and QED-HF files

For each QED-CCSD(2,2) or QED-HF text summary, a normalized CSV with the same
stem is available:

| File pattern | Contents |
| --- | --- |
| `relaxed_dir_<theta>_<phi>_QEDCCSD22.csv` | QED-CCSD(2,2) energies at QED-DFT-relaxed geometries. |
| `unrelaxed_dir_<theta>_<phi>_QEDCCSD22.csv` | QED-CCSD(2,2) energies at unrelaxed geometries. |
| `relaxed_dir_<theta>_<phi>_QEDHF.csv` | QED-HF energies at QED-DFT-relaxed geometries. |
| `unrelaxed_dir_<theta>_<phi>_QEDHF.csv` | QED-HF energies at unrelaxed geometries. |

The `.txt` files are the original reduced summaries. The `.csv` files are the
preferred machine-readable summaries because they include method, geometry,
direction, lambda-vector components, absolute energies, and relative energies.

## pQED files

The `unrelaxed_dir_*_pqed_49_*.csv` files are pQED summaries for unrelaxed
geometries. They duplicate the files in `../../pqed/`, which are the clearer
method-level copies used by `../../scripts/plot_isomer_energies.py`.

## Relative energy sign convention

Relative energies are reported as `E(A) - E(B)` in kcal/mol:

| Direction | Column |
| --- | --- |
| `(theta=70, phi=31)` | `dE_ortho_meta_kcal/mol` |
| `(theta=65, phi=78)` | `dE_para_meta_kcal/mol` |
