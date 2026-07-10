#!/usr/bin/env python3
"""
Per-isomer-pair CS vs non-CS comparison maps.

For each Nel_X_Nph_Y run folder produced by pf_isomer_scan.py (auto-discovered
below this script; must contain both the CS and non-CS total-energy and
decomposition files), this produces two figures per isomer pair
(Ortho-Meta and Para-Meta):

  1. <tag>_total_energy_diff_CS_vs_nonCS_{pair}.png
     2 panels: CS (left) vs non-CS (right) total-energy difference.

  2. <tag>_decomp_diff_CS_vs_nonCS_{pair}.png
     4 rows (electronic, photonic, bilinear coupling, dipole self-energy)
     x 2 columns (CS left, non-CS right) of the decomposed-energy difference.

Reuses the column layout / plotting conventions of deploy_plots.py.
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Configuration ----------------------------------------------------------
# Folders are auto-discovered: any Nel_*_Nph_* directory (relative to this
# script) that contains both the CS and non-CS total-energy/decomposition
# files produced by pf_isomer_scan.py.  Set FOLDERS to a non-empty list to
# restrict/override which folders get plotted.
FOLDERS = []

DPI = 300
AU_TO_KCAL = 627.509
CMAP = "RdBu_r"

TOTAL_COLS = {
    "theta": 0, "phi": 1,
    "para_E": 5, "ortho_E": 6, "meta_E": 7,
}

DECOMP_COLS = {
    "theta": 0, "phi": 1,
    "para_E_el": 5,  "para_E_ph": 6,  "para_E_blc": 7,  "para_E_dse": 8,
    "ortho_E_el": 9, "ortho_E_ph": 10, "ortho_E_blc": 11, "ortho_E_dse": 12,
    "meta_E_el": 13, "meta_E_ph": 14, "meta_E_blc": 15, "meta_E_dse": 16,
}

# component key -> panel title
COMPONENTS = [
    ("el",  "Electronic"),
    ("ph",  "Photonic"),
    ("blc", "Bilinear Coupling"),
    ("dse", "Dipole Self-Energy"),
]


def discover_folders():
    """Nel_*_Nph_* folders (relative to this script) with both CS and
    non-CS total-energy/decomposition files."""
    if FOLDERS:
        return FOLDERS
    found = []
    for path in sorted(glob.glob("Nel_*_Nph_*")):
        tag = os.path.basename(path)
        needed = [
            f"isomer_{tag}_total_energies.dat",
            f"isomer_{tag}_total_energies_CS.dat",
            f"isomer_{tag}_energy_decomposition.dat",
            f"isomer_{tag}_energy_decomposition_CS.dat",
        ]
        if all(os.path.exists(os.path.join(path, n)) for n in needed):
            found.append(path)
    return found


def _grid(data, cols):
    """theta/phi meshgrid-style arrays reshaped from the flat scan data."""
    theta = data[:, cols["theta"]]
    phi = data[:, cols["phi"]]
    num_t = len(np.unique(theta))
    num_p = len(np.unique(phi))
    T = theta.reshape(num_t, num_p)
    P = phi.reshape(num_t, num_p)
    return T, P, num_t, num_p


def _style_axis(ax, title):
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(r"Azimuthal Angle $\phi$ (deg)", fontsize=11)
    ax.set_ylabel(r"Polar Angle $\theta$ (deg)", fontsize=11)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 45, 90, 135, 180])
    ax.grid(True, linestyle="--", alpha=0.3)


def plot_total_energy_pair(folder, tag, iso_a, iso_b):
    """2-panel CS (left) vs non-CS (right) total-energy difference (iso_a - iso_b)."""
    noncs = np.genfromtxt(os.path.join(folder, f"isomer_{tag}_total_energies.dat"), skip_header=2)
    cs = np.genfromtxt(os.path.join(folder, f"isomer_{tag}_total_energies_CS.dat"), skip_header=2)

    T, P, num_t, num_p = _grid(noncs, TOTAL_COLS)

    diff_cs = (cs[:, TOTAL_COLS[f"{iso_a}_E"]] - cs[:, TOTAL_COLS[f"{iso_b}_E"]]).reshape(num_t, num_p) * AU_TO_KCAL
    diff_noncs = (noncs[:, TOTAL_COLS[f"{iso_a}_E"]] - noncs[:, TOTAL_COLS[f"{iso_b}_E"]]).reshape(num_t, num_p) * AU_TO_KCAL
    print(F"min diff_cs is {diff_cs.min()}")
    print(F"max diff_cs is {diff_cs.max()}")
    print(F"min diff_noncs is {diff_noncs.min()}")
    print(F"max diff_noncs is {diff_noncs.max()}")
    # shared color scale so the two panels are directly comparable
    vmax = max(np.abs(diff_cs).max(), np.abs(diff_noncs).max())
    vmin = -vmax

    fig, (ax_cs, ax_noncs) = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)

    for ax, diff, title in ((ax_cs, diff_cs, "CS"), (ax_noncs, diff_noncs, "non-CS")):
        im = ax.pcolormesh(P, T, diff, shading="gouraud", cmap=CMAP) #, vmin=vmin, vmax=vmax)
        _style_axis(ax, title)
        fig.colorbar(im, ax=ax, label=r"$\Delta E$ (kcal / mol)")

    pair_pretty = f"{iso_a.capitalize()} $-$ {iso_b.capitalize()}"
    plt.suptitle(f"Total Energy Difference: {pair_pretty}  [{tag}]", fontsize=16, fontweight="bold")

    out = os.path.join(folder, f"{tag}_total_energy_diff_CS_vs_nonCS_{iso_a}-{iso_b}.png")
    plt.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"    saved {os.path.basename(out)}")


def plot_decomp_pair(folder, tag, iso_a, iso_b):
    """4 (component) x 2 (CS | non-CS) grid of decomposed-energy differences (iso_a - iso_b)."""
    noncs = np.genfromtxt(os.path.join(folder, f"isomer_{tag}_energy_decomposition.dat"), skip_header=2)
    cs = np.genfromtxt(os.path.join(folder, f"isomer_{tag}_energy_decomposition_CS.dat"), skip_header=2)

    T, P, num_t, num_p = _grid(noncs, DECOMP_COLS)

    fig, axes = plt.subplots(4, 2, figsize=(13, 18), constrained_layout=True)

    for row, (comp, title) in enumerate(COMPONENTS):
        diff_cs = (cs[:, DECOMP_COLS[f"{iso_a}_E_{comp}"]] - cs[:, DECOMP_COLS[f"{iso_b}_E_{comp}"]]).reshape(num_t, num_p) * AU_TO_KCAL
        diff_noncs = (noncs[:, DECOMP_COLS[f"{iso_a}_E_{comp}"]] - noncs[:, DECOMP_COLS[f"{iso_b}_E_{comp}"]]).reshape(num_t, num_p) * AU_TO_KCAL

        # shared color scale within a row so CS vs non-CS are directly comparable
        vmax = max(np.abs(diff_cs).max(), np.abs(diff_noncs).max())
        vmin = -vmax

        ax_cs, ax_noncs = axes[row, 0], axes[row, 1]
        im_cs = ax_cs.pcolormesh(P, T, diff_cs, shading="gouraud", cmap=CMAP, vmin=vmin, vmax=vmax)
        im_noncs = ax_noncs.pcolormesh(P, T, diff_noncs, shading="gouraud", cmap=CMAP, vmin=vmin, vmax=vmax)

        _style_axis(ax_cs, f"{title} (CS)")
        _style_axis(ax_noncs, f"{title} (non-CS)")
        fig.colorbar(im_cs, ax=ax_cs, label=r"$\Delta E$ (kcal / mol)")
        fig.colorbar(im_noncs, ax=ax_noncs, label=r"$\Delta E$ (kcal / mol)")

    pair_pretty = f"{iso_a.capitalize()} $-$ {iso_b.capitalize()}"
    plt.suptitle(f"Energy-Component Decomposition: {pair_pretty}  [{tag}]", fontsize=16, fontweight="bold")

    out = os.path.join(folder, f"{tag}_decomp_diff_CS_vs_nonCS_{iso_a}-{iso_b}.png")
    plt.savefig(out, dpi=DPI)
    plt.close(fig)
    print(f"    saved {os.path.basename(out)}")


def process_folder(folder):
    tag = folder
    print(f"[plot] {tag}")
    for iso_a, iso_b in (("ortho", "meta"), ("para", "meta")):
        plot_total_energy_pair(folder, tag, iso_a, iso_b)
        plot_decomp_pair(folder, tag, iso_a, iso_b)


def main():
    folders = discover_folders()
    if not folders:
        print("No Nel_*_Nph_* folders with both CS and non-CS output files were found.")
        return
    for folder in folders:
        process_folder(folder)


if __name__ == "__main__":
    main()
