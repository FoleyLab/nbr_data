"""
Publication-style QED-CCSD relaxed/unrelaxed comparison plot.

This is a focused companion to plot_isomer_energies.py. It reads the normalized
QED-CCSD(2,2) CSV summaries in Lambda_Scan_Results/QED_CCSD/summary and creates one figure with
four curves:

  - ortho - meta, theta=70, phi=31, relaxed geometry
  - ortho - meta, theta=70, phi=31, unrelaxed geometry
  - para - meta, theta=65, phi=78, relaxed geometry
  - para - meta, theta=65, phi=78, unrelaxed geometry

Visual conventions are borrowed from the deprecated final_plot_script.py:
coupling strength on the lower x axis, mode volume on the upper x axis, and a
thermal stabilization guide line. Here the guide is -5 k_B T at 298 K rather
than -10 k_B T, because the new scan is less strongly stabilizing.

Run:
    python plot_qed_ccsd_relaxed_unrelaxed_styled.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARY_DIR = SCRIPT_DIR / "Lambda_Scan_Results" / "QED_CCSD" / "summary"

OUTPUT_STEM = "qed_ccsd_relaxed_vs_unrelaxed_styled"

# Mode-volume conversion used in DEPRICATED_LAMBDA_SCAN_DATA/final_plot_script.py.
A0_NM = 0.05291772109
C_NM3 = 4.0 * math.pi * (A0_NM**3)

# Thermal reference: R*T in kcal/mol at 298.15 K.
RT_KCAL_MOL_298K = 0.00198720425864083 * 298.15
THERMAL_THRESHOLD_KBT_A = 2
THERMAL_THRESHOLD_KBT_B = 4
THERMAL_THRESHOLD_KBT_C = 6

THERMAL_THRESHOLD_KCAL_A = -THERMAL_THRESHOLD_KBT_A * RT_KCAL_MOL_298K
THERMAL_THRESHOLD_KCAL_B = -THERMAL_THRESHOLD_KBT_B * RT_KCAL_MOL_298K
THERMAL_THRESHOLD_KCAL_C = -THERMAL_THRESHOLD_KBT_C * RT_KCAL_MOL_298K

# Optional one-sided model-bias envelope from cavity-free CCSD(T)-CCSD shifts.
# This is not a statistical error bar.
SHOW_MISSING_TRIPLES_ENVELOPE = False
MISSING_TRIPLES_STABILIZATION_KCAL = {
    "ortho_meta": 1.0,
    "para_meta": 1.5,
}
MISSING_TRIPLES_ALPHA = 0.13

COLORS = {
    "ortho_meta": "#0072B2",  # blue, color-blind friendly
    "para_meta": "#D55E00",   # vermillion, color-blind friendly
}

SERIES = [
    {
        "geometry": "relaxed",
        "pair": "ortho_meta",
        "theta": 70,
        "phi": 31,
        "marker": "o",
        "label": (
            r"$\Delta E_{\mathrm{ortho-meta}},\; \theta=70^{\circ},\; "
            r"\phi=31^{\circ}$, relaxed"
        ),
    },
    {
        "geometry": "unrelaxed",
        "pair": "ortho_meta",
        "theta": 70,
        "phi": 31,
        "marker": "o",
        "label": (
            r"$\Delta E_{\mathrm{ortho-meta}},\; \theta=70^{\circ},\; "
            r"\phi=31^{\circ}$, unrelaxed"
        ),
    },
    {
        "geometry": "relaxed",
        "pair": "para_meta",
        "theta": 65,
        "phi": 78,
        "marker": "s",
        "label": (
            r"$\Delta E_{\mathrm{para-meta}},\; \theta=65^{\circ},\; "
            r"\phi=78^{\circ}$, relaxed"
        ),
    },
    {
        "geometry": "unrelaxed",
        "pair": "para_meta",
        "theta": 65,
        "phi": 78,
        "marker": "s",
        "label": (
            r"$\Delta E_{\mathrm{para-meta}},\; \theta=65^{\circ},\; "
            r"\phi=78^{\circ}$, unrelaxed"
        ),
    },
]


def lambda_to_mode_volume(lambda_au: float) -> float:
    """Convert coupling strength lambda in a.u. to mode volume in nm^3.

    Matplotlib's secondary axis machinery calls this with both scalars and
    NumPy arrays, so the implementation must be vector-safe.
    """
    lambda_au = np.asarray(lambda_au)
    out = np.full_like(lambda_au, np.inf, dtype=float)
    np.divide(C_NM3, lambda_au**2, out=out, where=lambda_au != 0)
    return out


def mode_volume_to_lambda(volume_nm3: float) -> float:
    """Inverse transform for matplotlib's secondary_xaxis."""
    volume_nm3 = np.asarray(volume_nm3)
    ratio = np.full_like(volume_nm3, np.inf, dtype=float)
    np.divide(C_NM3, volume_nm3, out=ratio, where=volume_nm3 != 0)
    return np.sqrt(ratio)


def read_qed_ccsd_series(entry: dict) -> tuple[list[float], list[float]]:
    """Read one normalized QED-CCSD summary CSV into sorted x/y arrays."""
    path = SUMMARY_DIR / (
        f"{entry['geometry']}_dir_{entry['theta']}_{entry['phi']}_QEDCCSD22.csv"
    )
    delta_column = f"dE_{entry['pair']}_kcal/mol"

    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    rows.sort(key=lambda row: float(row["lambda_magnitude"]))
    lambdas = [float(row["lambda_magnitude"]) for row in rows]
    delta_e = [float(row[delta_column]) for row in rows]
    return lambdas, delta_e


def configure_style() -> None:
    """Set plotting defaults once for a clean publication-style figure."""
    plt.rcParams.update(
        {
            "font.size": 22,
            "font.family": "sans-serif",
            "axes.labelsize": 22,
            "axes.titlesize": 22,
            "axes.linewidth": 1.2,
            "xtick.labelsize": 22,
            "ytick.labelsize": 22,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": True,
            "legend.fontsize": 18,
            "savefig.bbox": "tight",
        }
    )


def plot_series(ax) -> None:
    """Plot all configured QED-CCSD relaxed/unrelaxed curves."""
    for entry in SERIES:
        lambdas, delta_e = read_qed_ccsd_series(entry)
        color = COLORS[entry["pair"]]
        is_relaxed = entry["geometry"] == "relaxed"

        if SHOW_MISSING_TRIPLES_ENVELOPE:
            shift = MISSING_TRIPLES_STABILIZATION_KCAL[entry["pair"]]
            lambda_array = np.asarray(lambdas)
            delta_array = np.asarray(delta_e)
            ax.fill_between(
                lambda_array,
                delta_array - shift,
                delta_array,
                color=color,
                alpha=MISSING_TRIPLES_ALPHA,
                linewidth=0,
                zorder=1,
            )

        ax.plot(
            lambdas,
            delta_e,
            color=color,
            linestyle="-" if is_relaxed else "--",
            linewidth=2.6,
            marker=entry["marker"],
            markersize=8,
            markerfacecolor=color if is_relaxed else "white",
            markeredgecolor=color,
            markeredgewidth=1.8,
            label=entry["label"],
            zorder=2,
        )


def add_mode_volume_axis(ax) -> None:
    """Add mode volume as an upper x axis aligned to the lambda ticks."""
    secax = ax.secondary_xaxis(
        "top",
        functions=(lambda_to_mode_volume, mode_volume_to_lambda),
    )
    secax.set_xlabel(r"Mode Volume $V$ (nm$^3$)", labelpad=15)

    lambda_ticks = [0.02, 0.04, 0.06, 0.08, 0.10]
    volume_ticks = [lambda_to_mode_volume(value) for value in lambda_ticks]
    secax.set_ticks(volume_ticks)
    secax.set_xticklabels([f"{value:.2f}" for value in volume_ticks])


def add_thermal_threshold(ax) -> None:
    """Draw and label the -X k_B T stabilization guide line."""
    ax.axhline(
        THERMAL_THRESHOLD_KCAL_A,
        color="black",
        linestyle=":",
        linewidth=1.8,
        alpha=0.75,
        zorder=0,
    )
    ax.text(
        0.021,
        THERMAL_THRESHOLD_KCAL_A + 0.18,
        rf"$-{THERMAL_THRESHOLD_KBT_A}\;k_B T$ stabilization",
        fontsize=16,
        fontweight="bold",
        color="#333333",
    )
    ax.axhline(
        THERMAL_THRESHOLD_KCAL_B,
        color="black",
        linestyle=":",
        linewidth=1.8,
        alpha=0.75,
        zorder=0,
    )
    ax.text( 
        0.021,
        THERMAL_THRESHOLD_KCAL_B + 0.18,
        rf"$-{THERMAL_THRESHOLD_KBT_B}\;k_B T$ stabilization",
        fontsize=16,
        fontweight="bold",
        color="#333333",
    ) 
    ax.axhline(
        THERMAL_THRESHOLD_KCAL_C,
        color="black",
        linestyle=":",
        linewidth=1.8,
        alpha=0.75,
        zorder=0,
    )
    ax.text(
        0.021,
        THERMAL_THRESHOLD_KCAL_C + 0.18,
        rf"$-{THERMAL_THRESHOLD_KBT_C}\;k_B T$ stabilization",
        fontsize=16,
        fontweight="bold",
        color="#333333",
    )



def main() -> None:
    configure_style()
    fig, ax = plt.subplots(figsize=(11, 7.5), dpi=300)

    plot_series(ax)
    add_thermal_threshold(ax)
    add_mode_volume_axis(ax)

    ax.axhline(0.0, color="0.45", linestyle="-", linewidth=1.0, alpha=0.5, zorder=0)
    ax.set_xlabel(r"Coupling Strength $\lambda$ (a.u.)", labelpad=10)
    ax.set_ylabel(r"$\Delta E$ (kcal/mol)", labelpad=10)
    ax.set_xlim(0.015, 0.105)
    ax.set_ylim(-5.2, 5.4)
    ax.set_xticks([0.02, 0.04, 0.06, 0.08, 0.10])
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    #ax.legend(loc="lower left", frameon=True, shadow=True)
    handles, labels = ax.get_legend_handles_labels()
    if SHOW_MISSING_TRIPLES_ENVELOPE:
        handles.append(
            Patch(facecolor="0.35", alpha=MISSING_TRIPLES_ALPHA, edgecolor="none")
        )
        labels.append(r"estimated missing-(T) stabilization envelope")
    #ax.legend(handles, labels, bbox_to_anchor=(1.05,1), loc="upper left", frameon=True, shadow=True, borderaxespad=0)
    # Replaced ax.legend(...) call:
    ax.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.22),  # Positions legend above the top x-axis
        ncol=2,                       # Split into 2 columns to save vertical space
        frameon=True,
        shadow=False,
        fontsize=16,                  # Scale legend font slightly relative to 20pt axes
        borderaxespad=0,
    )

    # Use constrained_layout instead of tight_layout to handle secondary axes & legends gracefully
    fig.set_layout_engine("constrained")
    fig.tight_layout()
    png_path = SCRIPT_DIR / f"{OUTPUT_STEM}.png"
    pdf_path = SCRIPT_DIR / f"{OUTPUT_STEM}.pdf"
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Saved {png_path.name} and {pdf_path.name}")


if __name__ == "__main__":
    main()
