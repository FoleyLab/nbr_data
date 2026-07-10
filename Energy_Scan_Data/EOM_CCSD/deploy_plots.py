#!/usr/bin/env python3
"""
Deploy the multipanel difference plots across every Nel_X_Nph_Y folder.

For each folder that contains the three isomer .dat files, this script writes
PNGs directly into that folder:

  Decomposition maps (from ..._energy_decomposition.dat):
    <tag>_decomp_diff_el.png     electronic component  (Ortho-Meta, Para-Meta)
    <tag>_decomp_diff_ph.png     photonic component
    <tag>_decomp_diff_blc.png    bilinear-coupling component
    <tag>_decomp_diff_dse.png    dipole-self-energy component
    <tag>_decomp_diff_total_sum_of_components.png
                                 sum of the four component expectation values

  Total-energy map (from ..._total_energies.dat):
    <tag>_total_diff_lowest_eigenvalue.png
                                 raw ground-state total (lowest PF eigenvalue)

where <tag> is the folder name (e.g. Nel_10_Nph_10).

The two "total" files are named differently on purpose:
  *_total_sum_of_components*  -> sum of <H_el>+<H_ph>+<H_blc>+<H_dse>
  *_lowest_eigenvalue*        -> the diagonalized ground-state energy
Comparing them verifies the lowest eigenvalue equals the sum of the component
expectation values.  This script also prints the numeric residual between the
two for each folder (should be ~0 to machine precision).

Refactored from plot_multipanel_energy_decomposition.py and
plot_multipanel_ortho_para_meta.py.
"""

import os
import re
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless: save PNGs, never open a window
import matplotlib.pyplot as plt

# --- Configuration ----------------------------------------------------------
DPI = 300
AU_TO_KCAL = 627.509
CMAP = "RdBu_r"

# Column layout of the decomposition .dat files
DECOMP_COLS = {
    "theta": 0, "phi": 1,
    "para_E_el": 5,  "para_E_ph": 6,  "para_E_blc": 7,  "para_E_dse": 8,
    "ortho_E_el": 9, "ortho_E_ph": 10, "ortho_E_blc": 11, "ortho_E_dse": 12,
    "meta_E_el": 13, "meta_E_ph": 14, "meta_E_blc": 15, "meta_E_dse": 16,
}

# Column layout of the total-energy .dat files
TOTAL_COLS = {"theta": 0, "phi": 1, "para": 5, "ortho": 6, "meta": 7}


def _make_panel(theta, phi, e_para, e_ortho, e_meta,
                legend_label, suptitle, out_path):
    """Two-panel (Ortho-Meta, Para-Meta) pcolormesh of an energy difference."""
    num_t = len(np.unique(theta))
    num_p = len(np.unique(phi))

    T = theta.reshape(num_t, num_p)
    P = phi.reshape(num_t, num_p)

    diff_ortho_meta = (e_ortho - e_meta).reshape(num_t, num_p) * AU_TO_KCAL
    diff_para_meta = (e_para - e_meta).reshape(num_t, num_p) * AU_TO_KCAL

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    im1 = ax1.pcolormesh(P, T, diff_ortho_meta, shading="gouraud", cmap=CMAP)
    ax1.set_title(r"$\Delta E$ (Ortho $-$ Meta)", fontsize=14)
    ax1.set_xlabel(r"Azimuthal Angle $\phi$ (deg)", fontsize=12)
    ax1.set_ylabel(r"Polar Angle $\theta$ (deg)", fontsize=12)
    fig.colorbar(im1, ax=ax1, label=legend_label)

    im2 = ax2.pcolormesh(P, T, diff_para_meta, shading="gouraud", cmap=CMAP)
    ax2.set_title(r"$\Delta E$ (Para $-$ Meta)", fontsize=14)
    ax2.set_xlabel(r"Azimuthal Angle $\phi$ (deg)", fontsize=12)
    ax2.set_ylabel(r"Polar Angle $\theta$ (deg)", fontsize=12)
    fig.colorbar(im2, ax=ax2, label=legend_label)

    for ax in (ax1, ax2):
        ax.set_xticks([0, 90, 180, 270, 360])
        ax.set_yticks([0, 45, 90, 135, 180])
        ax.grid(True, linestyle="--", alpha=0.3)

    plt.suptitle(suptitle, fontsize=16, fontweight="bold")
    plt.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f"    saved {os.path.basename(out_path)}")


def plot_decomposition(decomp_file, tag, folder):
    """Generate the per-component and sum-of-components difference maps."""
    data = np.genfromtxt(decomp_file, skip_header=2)
    theta = data[:, DECOMP_COLS["theta"]]
    phi = data[:, DECOMP_COLS["phi"]]

    components = {
        "el":  "Electronic Energy Difference (kcal / mol)",
        "ph":  "Photonic Energy Difference (kcal / mol)",
        "blc": "Bilinear Coupling Energy Difference (kcal / mol)",
        "dse": "Dipole Self Energy Difference (kcal / mol)",
    }

    for comp, label in components.items():
        e_para = data[:, DECOMP_COLS[f"para_E_{comp}"]]
        e_ortho = data[:, DECOMP_COLS[f"ortho_E_{comp}"]]
        e_meta = data[:, DECOMP_COLS[f"meta_E_{comp}"]]
        out = os.path.join(folder, f"{tag}_decomp_diff_{comp}.png")
        _make_panel(theta, phi, e_para, e_ortho, e_meta,
                    label, f"{label}  [{tag}]", out)

    # sum of the four component expectation values
    def total(iso):
        return (data[:, DECOMP_COLS[f"{iso}_E_el"]]
                + data[:, DECOMP_COLS[f"{iso}_E_ph"]]
                + data[:, DECOMP_COLS[f"{iso}_E_blc"]]
                + data[:, DECOMP_COLS[f"{iso}_E_dse"]])

    sum_para, sum_ortho, sum_meta = total("para"), total("ortho"), total("meta")
    label = "Total Energy Difference (Sum of Components, kcal / mol)"
    out = os.path.join(folder, f"{tag}_decomp_diff_total_sum_of_components.png")
    _make_panel(theta, phi, sum_para, sum_ortho, sum_meta,
                label, f"Sum-of-Components Total  [{tag}]", out)

    return sum_para, sum_ortho, sum_meta


def plot_total(total_file, tag, folder):
    """Generate the raw total-energy (lowest eigenvalue) difference map."""
    data = np.genfromtxt(total_file, skip_header=2)
    theta = data[:, TOTAL_COLS["theta"]]
    phi = data[:, TOTAL_COLS["phi"]]
    e_para = data[:, TOTAL_COLS["para"]]
    e_ortho = data[:, TOTAL_COLS["ortho"]]
    e_meta = data[:, TOTAL_COLS["meta"]]

    label = "Total Energy Difference (Lowest Eigenvalue, kcal / mol)"
    out = os.path.join(folder, f"{tag}_total_diff_lowest_eigenvalue.png")
    _make_panel(theta, phi, e_para, e_ortho, e_meta,
                label, f"Lowest-Eigenvalue Total  [{tag}]", out)

    return e_para, e_ortho, e_meta


def find_files(folder, tag):
    """Locate the decomposition and total .dat files inside a folder."""
    decomp = os.path.join(folder, f"isomer_{tag}_energy_decomposition.dat")
    total = os.path.join(folder, f"isomer_{tag}_total_energies.dat")
    # fall back to any matching file if the exact tag name differs
    if not os.path.exists(decomp):
        hits = glob.glob(os.path.join(folder, "*decomposition*.dat"))
        decomp = hits[0] if hits else None
    if not os.path.exists(total):
        hits = glob.glob(os.path.join(folder, "*total_energies*.dat"))
        total = hits[0] if hits else None
    return decomp, total


def main():
    folders = sorted(
        d for d in glob.glob("Nel_*_Nph_*")
        if os.path.isdir(d) and re.match(r"Nel_\d+_Nph_\d+$", d)
    )
    if not folders:
        print("No Nel_X_Nph_Y folders found.")
        return

    print(f"Found {len(folders)} candidate folders.\n")
    residual_report = []

    for folder in folders:
        tag = folder.rstrip("/")
        decomp, total = find_files(folder, tag)

        if not decomp or not total:
            print(f"[skip] {tag}: missing .dat files")
            continue

        print(f"[plot] {tag}")
        sum_p, sum_o, sum_m = plot_decomposition(decomp, tag, folder)
        raw_p, raw_o, raw_m = plot_total(total, tag, folder)

        # verify lowest eigenvalue == sum of component expectation values
        res = max(
            np.max(np.abs(raw_p - sum_p)),
            np.max(np.abs(raw_o - sum_o)),
            np.max(np.abs(raw_m - sum_m)),
        )
        residual_report.append((tag, res))
        print(f"    max|eigenvalue - sum_of_components| = {res:.3e} Ha\n")

    # summary
    print("=" * 60)
    print("Residual check (lowest eigenvalue vs. sum of components):")
    print(f"{'folder':<18}{'max abs residual (Ha)':>24}")
    print("-" * 42)
    for tag, res in residual_report:
        flag = "" if res < 1e-9 else "   <-- NONZERO"
        print(f"{tag:<18}{res:>24.3e}{flag}")


if __name__ == "__main__":
    main()
