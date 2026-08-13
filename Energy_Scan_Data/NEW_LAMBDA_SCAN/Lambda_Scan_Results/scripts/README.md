# scripts - Analysis and plotting tools

| File | Purpose |
| --- | --- |
| `plot_isomer_energies.py` | Generates method-comparison plots of relative isomer energies vs. `|lambda|`. |
| `plot_qed_ccsd_relaxed_unrelaxed_styled.py` | Focused publication-style QED-CCSD relaxed/unrelaxed comparison with lower coupling-strength axis, upper mode-volume axis, and `-5 k_B T` guide line. |

## plot_isomer_energies.py

Run from this folder or from anywhere:

```bash
python plot_isomer_energies.py
```

The script writes both `.png` and `.pdf` outputs into `scripts/`. The default
configuration regenerates these comparison figures:

| Figure stem | Comparison |
| --- | --- |
| `01_relaxed_vs_unrelaxed_no_zpe` | QED-DFT relaxed vs. unrelaxed for both target pairs. |
| `02_relaxed_zpe_vs_no_zpe` | Relaxed QED-DFT raw electronic vs. ZPE-corrected relative energies. |
| `03_qeddft_vs_pqed_ortho_meta` | pQED vs. QED-DFT relaxed and unrelaxed for ortho-meta at `(70,31)`. |
| `04_qeddft_vs_pqed_para_meta` | pQED vs. QED-DFT relaxed and unrelaxed for para-meta at `(65,78)`. |
| `05_qed_ccsd_relaxed_vs_unrelaxed` | QED-CCSD relaxed vs. unrelaxed for both target pairs. |
| `06_pqed_vs_qed_ccsd_ortho_meta` | pQED vs. QED-CCSD relaxed and unrelaxed for ortho-meta at `(70,31)`. |
| `07_pqed_vs_qed_ccsd_para_meta` | pQED vs. QED-CCSD relaxed and unrelaxed for para-meta at `(65,78)`. |
| `08_pqed_coherent_state_check` | pQED coherent-state vs. non-coherent-state comparison. |

## Data read by the script

| Method key | Reads from | Notes |
| --- | --- | --- |
| `qed_ccsd` | `../QED_CCSD/summary/*_QEDCCSD22.csv` | QED-CCSD(2,2), relaxed or unrelaxed. |
| `qed_dft` | `../qed_dft_relaxed/*.csv`, `../qed_dft_unrelaxed/*.csv` | QED-DFT, relaxed or unrelaxed; relaxed data can use `zpe=True`. |
| `pqed` | `../pqed/pqed_49_*_scan*.csv` | pQED unrelaxed data; supports `Nph=3/10` and `CS=True/False`. |

To add a plot, edit the `DEFAULT_FIGURES` dictionary near the bottom of the
script. Each curve is a small dictionary describing `method`, `pair`,
`direction`, and any method-specific options.

Requires Python with `matplotlib`. CSV parsing uses the Python standard
library.

## plot_qed_ccsd_relaxed_unrelaxed_styled.py

Run:

```bash
python plot_qed_ccsd_relaxed_unrelaxed_styled.py
```

This script reads only the normalized QED-CCSD(2,2) summaries in
`../QED_CCSD/summary/` and writes:

| Output | Contents |
| --- | --- |
| `qed_ccsd_relaxed_vs_unrelaxed_styled.png` | 300 dpi raster figure. |
| `qed_ccsd_relaxed_vs_unrelaxed_styled.pdf` | Vector figure. |

The plot uses blue for `ortho - meta` at `(theta=70, phi=31)` and vermillion
for `para - meta` at `(theta=65, phi=78)`. Relaxed geometries are solid with
filled markers; unrelaxed geometries are dashed with open markers.
