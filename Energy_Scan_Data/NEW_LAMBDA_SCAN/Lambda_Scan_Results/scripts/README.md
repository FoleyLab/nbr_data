# scripts — Analysis and plotting tools

| File | Purpose |
| --- | --- |
| `plot_isomer_energies.py` | Publication-quality plots of isomer-pair energy differences (kcal/mol) vs. \|λ\| for the QED-DFT / pQED data in the sibling folders. Config-driven: edit the `FIGURES` dict, run `python plot_isomer_energies.py`, and it emits `.png` (300 dpi) + `.pdf` figures in this folder. |

## plot_isomer_energies.py

Reads data from the sibling folders of this script:

| Method | Reads from |
| --- | --- |
| `qed_dft_relaxed` | `../qed_dft_relaxed/*_hartree.csv` (also used for the +ZPE curves) |
| `qed_dft_no_relax` | `../qed_dft_unrelaxed/*_scan_qed_dft_no_relax.csv` |
| `pqed` | `../pqed/pqed_49_*_scan[_CS].csv` |

See the module docstring for the exact `FIGURES` entry syntax. Running
`python plot_isomer_energies.py` regenerates the current four figures:

| Figure | Comparison |
| --- | --- |
| `01_relaxed_vs_unrelaxed_no_zpe` | relaxed vs unrelaxed QED-DFT |
| `02_relaxed_zpe_vs_no_zpe` | relaxed QED-DFT with vs without ZPE |
| `03_qeddft_vs_pqed_ortho_meta` | QED-DFT vs pQED, ortho−meta |
| `04_qeddft_vs_pqed_para_meta` | QED-DFT vs pQED, para−meta |

Figures (`.png` 300 dpi + `.pdf`) are saved in this folder (`scripts/`).
Requires `python` with `numpy`, `pandas`, `matplotlib` (e.g. the `jbook`
conda env).
