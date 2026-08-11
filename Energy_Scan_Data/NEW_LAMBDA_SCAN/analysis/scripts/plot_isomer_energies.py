"""
plot_isomer_energies.py

Publication-quality plots of isomer-pair energy differences (kcal/mol) vs.
cavity-coupling magnitude |lambda|, for the nitrobenzene QED-DFT / pQED
scan data in this analysis folder.

DATA LAYOUT (sibling folders of this script):
  qed_dft_relaxed/    -- relaxed QED-DFT per-lambda tables (Hartree), from analyze.py
  qed_dft_unrelaxed/  -- unrelaxed QED-DFT scans (kcal/mol)
  pqed/               -- parameterized QED scans (kcal/mol)

============================== HOW TO USE ==================================
Everything you normally need to change lives in the FIGURES dict near the
bottom of this file. Each entry is:

    "output_file_stem": [ series_1, series_2, ... ]

and each series is a small dict describing ONE curve. Put more than one
series dict in a figure's list to overlay multiple methods/directions on
the same axes. Three "method" kinds are supported:

  qed_dft_relaxed   -- fully relaxed geometry, QED-DFT
      requires: pair ("ortho_meta" or "para_meta"), direction (theta, phi)
      optional: zpe (True/False, default False) -- include ZPE correction

  qed_dft_no_relax  -- unrelaxed geometry, QED-DFT
      requires: pair ("ortho_meta" or "para_meta"), direction (theta, phi)
      (no ZPE is available/meaningful for unrelaxed geometries)

  pqed              -- unrelaxed geometry, parameterized QED
      requires: pair ("ortho_meta", "para_meta", or "ortho_para"),
                direction (theta, phi), Nel, Nph
      optional: CS (True/False, default False) -- coherent-state transform

Every series dict may also include an optional "label" string to override
the auto-generated legend text.

Run this file directly (`python plot_isomer_energies.py`) to regenerate
every figure in FIGURES as both .png (300 dpi) and .pdf in this folder.
==============================================================================
"""

from pathlib import Path
from itertools import cycle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------
# Paths / constants
# -----------------------------------------------------------------------

DATA_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = DATA_DIR  # where generated figures are saved

RELAXED_DIR = DATA_DIR.parent / "qed_dft_relaxed"
UNRELAXED_DIR = DATA_DIR.parent / "qed_dft_unrelaxed"
PQED_DIR = DATA_DIR.parent / "pqed"

HARTREE_TO_KCAL = 627.5094740631

PAIR_LABELS = {
    "ortho_meta": "ortho − meta",
    "para_meta": "para − meta",
    "ortho_para": "ortho − para",
}

# -----------------------------------------------------------------------
# Publication-style rcParams (tweak once, applies everywhere)
# -----------------------------------------------------------------------

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "axes.linewidth": 1.0,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "lines.linewidth": 1.8,
    "lines.markersize": 6,
    "savefig.bbox": "tight",
})

MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*"]


# -----------------------------------------------------------------------
# Loaders -- one per data "kind". Each returns (lambda_array, dE_kcal_array)
# -----------------------------------------------------------------------

def load_qed_dft_relaxed(pair, direction, zpe=False):
    theta, phi = direction
    if pair == "ortho_meta":
        fname = f"ortho-meta_{theta}_{phi}_hartree.csv"
    elif pair == "para_meta":
        fname = f"para-meta_{theta}_{phi}_hartree.csv"
    else:
        raise ValueError(
            f"qed_dft_relaxed only supports pair='ortho_meta' or 'para_meta', got {pair!r}"
        )

    path = RELAXED_DIR / fname
    if not path.exists():
        raise FileNotFoundError(f"Expected relaxed QED-DFT file not found: {path.name}")

    df = pd.read_csv(path, comment="#")
    lam = df["magnitude"].to_numpy(dtype=float)
    col = "dE_zpe_Ha" if zpe else "dE_raw_Ha"
    dE = df[col].to_numpy(dtype=float) * HARTREE_TO_KCAL
    return lam, dE


def load_qed_dft_no_relax(pair, direction):
    theta, phi = direction
    if pair == "ortho_meta":
        fname = f"ortho_meta_dir{theta}_{phi}_scan_qed_dft_no_relax.csv"
    elif pair == "para_meta":
        fname = f"para_meta_dir{theta}_{phi}_scan_qed_dft_no_relax.csv"
    else:
        raise ValueError(
            f"qed_dft_no_relax only supports pair='ortho_meta' or 'para_meta', got {pair!r}"
        )

    path = UNRELAXED_DIR / fname
    if not path.exists():
        raise FileNotFoundError(f"Expected unrelaxed QED-DFT file not found: {path.name}")

    df = pd.read_csv(path, comment="#")
    lam = np.sqrt(df["Ex"] ** 2 + df["Ey"] ** 2 + df["Ez"] ** 2).to_numpy(dtype=float)
    dE = df["dE_kcal/mol"].to_numpy(dtype=float)
    return lam, dE


def load_pqed(pair, direction, Nel, Nph, CS=False):
    theta, phi = direction
    if pair not in PAIR_LABELS:
        raise ValueError(f"Unknown pair {pair!r}")

    cs_suffix = "_CS" if CS else ""
    fname = f"pqed_{Nel}_{Nph}_dir{theta}_{phi}_scan{cs_suffix}.csv"
    path = PQED_DIR / fname
    if not path.exists():
        raise FileNotFoundError(f"Expected pQED file not found: {path.name}")

    df = pd.read_csv(path, comment="#")
    lam = df["lambda_magnitude"].to_numpy(dtype=float)
    dE = df[f"dE_{pair}_kcal/mol"].to_numpy(dtype=float)
    return lam, dE


# -----------------------------------------------------------------------
# Dispatch + auto-labeling
# -----------------------------------------------------------------------

def _auto_label(entry):
    theta, phi = entry["direction"]
    pair_txt = PAIR_LABELS[entry["pair"]]
    dir_txt = f"θ={theta}° ϕ={phi}°"
    method = entry["method"]

    if method == "qed_dft_relaxed":
        zpe_txt = "+ZPE" if entry.get("zpe", False) else "no ZPE"
        return f"QED-DFT relaxed, {zpe_txt} ({pair_txt}), {dir_txt}"
    elif method == "qed_dft_no_relax":
        return f"QED-DFT unrelaxed ({pair_txt}), {dir_txt}"
    elif method == "pqed":
        cs_txt = "CS" if entry.get("CS", False) else "no CS"
        return f"pQED ({entry['Nel']}e,{entry['Nph']}ph) {cs_txt} ({pair_txt}), {dir_txt}"
    else:
        raise ValueError(f"Unknown method {method!r}")


def get_series(entry):
    """Load one series dict -> (lambda, dE_kcal, legend_label)."""
    method = entry["method"]
    pair = entry["pair"]
    direction = entry["direction"]

    if method == "qed_dft_relaxed":
        lam, dE = load_qed_dft_relaxed(pair, direction, entry.get("zpe", False))
    elif method == "qed_dft_no_relax":
        lam, dE = load_qed_dft_no_relax(pair, direction)
    elif method == "pqed":
        lam, dE = load_pqed(pair, direction, entry["Nel"], entry["Nph"], entry.get("CS", False))
    else:
        raise ValueError(f"Unknown method {method!r}")

    label = entry.get("label") or _auto_label(entry)

    order = np.argsort(lam)
    return lam[order], dE[order], label


# -----------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------

def make_plot(series_list, output_name, title=None):
    """Overlay every series dict in series_list on one figure and save it."""
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    markers = cycle(MARKER_CYCLE)

    for entry in series_list:
        lam, dE, label = get_series(entry)
        ax.plot(lam, dE, marker=next(markers), label=label)

    ax.axhline(0.0, color="0.6", linewidth=0.9, linestyle="--", zorder=0)
    ax.set_xlabel(r"$|\lambda|$ (a.u.)")
    ax.set_ylabel(r"$\Delta E$ (kcal/mol)")
    if title:
        ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()

    png_path = OUTPUT_DIR / f"{output_name}.png"
    pdf_path = OUTPUT_DIR / f"{output_name}.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Saved {png_path.name} and {pdf_path.name}")


# =========================================================================
# USER CONFIGURATION -- edit this dict to choose what gets plotted.
#
# Each key is the output filename stem; each value is a list of series
# dicts (see the module docstring above for the fields each "method"
# needs). One series per list = a single-curve plot; multiple series =
# an overlay comparison plot.
#
# NOTE: direction (65, 78) only has para-meta data on disk (no ortho-meta
# relaxed/unrelaxed files exist for that direction), so the (65,78)
# examples below use "para_meta" rather than "ortho_meta".
# =========================================================================

FIGURES = {
    # 1. Relaxed vs Unrelaxed
    "01_relaxed_vs_unrelaxed_no_zpe": [
        dict(method="qed_dft_relaxed", pair="ortho_meta", direction=(70, 31), zpe=False),
        dict(method="qed_dft_no_relax", pair="ortho_meta", direction=(70, 31)),
        dict(method="qed_dft_relaxed", pair="para_meta", direction=(65, 78), zpe=False),
	dict(method="qed_dft_no_relax", pair="para_meta", direction=(65, 78)),
    ],

    # 2. Relaxed with and without ZPE
    "02_relaxed_zpe_vs_no_zpe": [
        dict(method="qed_dft_relaxed", pair="ortho_meta", direction=(70, 31), zpe=False),
        dict(method="qed_dft_relaxed", pair="ortho_meta", direction=(70, 31), zpe=True),
        dict(method="qed_dft_relaxed", pair="para_meta", direction=(65, 78), zpe=False),
        dict(method="qed_dft_relaxed", pair="para_meta", direction=(65, 78), zpe=True)

    ],
    "03_qeddft_vs_pqed_ortho_meta": [
        dict(method="qed_dft_relaxed", pair="ortho_meta", direction=(70, 31), zpe=True),
        dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=10, CS=False),
        dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=3, CS=False),

    ],
    "04_qeddft_vs_pqed_para_meta": [
        dict(method="qed_dft_relaxed", pair="para_meta", direction=(65, 78), zpe=True),
        dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=10, CS=False),
        dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=3, CS=False),
    ],
    # 3. Unrelaxed ortho-meta, QED-DFT, dir (70,31)
    #"03_ortho_meta_dir70_31_unrelaxed_qeddft": [
    #    dict(method="qed_dft_no_relax", pair="ortho_meta", direction=(70, 31)),
    #],

    # 4. Unrelaxed ortho-meta, pQED (Nel=49, Nph=10), with vs. without CS, dir (70,31)
    #"04_ortho_meta_dir70_31_pqed_CS_compare": [
    #    dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=10, CS=False),
    #    dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=10, CS=True),
    #],

    # 5. Relaxed para-meta, no ZPE, QED-DFT, dir (65,78)
    #"05_para_meta_dir65_78_relaxed_noZPE": [
    #    dict(method="qed_dft_relaxed", pair="para_meta", direction=(65, 78), zpe=False),
    #],

    # 6. Relaxed para-meta, with ZPE, QED-DFT, dir (65,78)
    #"06_para_meta_dir65_78_relaxed_ZPE": [
    #    dict(method="qed_dft_relaxed", pair="para_meta", direction=(65, 78), zpe=True),
    #],

    # 7. Unrelaxed para-meta, QED-DFT, dir (65,78)
    #"07_para_meta_dir65_78_unrelaxed_qeddft": [
    #    dict(method="qed_dft_no_relax", pair="para_meta", direction=(65, 78)),
    #],

    # 8. Unrelaxed para-meta, pQED (Nel=49, Nph=10), with vs. without CS, dir (65,78)
    #"08_para_meta_dir65_78_pqed_CS_compare": [
    #    dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=10, CS=False),
    #    dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=10, CS=True),
    #],
}


if __name__ == "__main__":
    for output_name, series_list in FIGURES.items():
        make_plot(series_list, output_name)
