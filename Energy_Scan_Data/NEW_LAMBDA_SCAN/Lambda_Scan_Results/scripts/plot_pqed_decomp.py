"""
Publication-style pQED energy-decomposition plots.

This script reads pQED decomposition CSVs and plots relative energy
differences by Hamiltonian component:

    Delta E_component = (E_A_component - E_B_component) * 627.509

For this dataset, the two default panels are:

  - ortho - meta at theta=70, phi=31
  - para - meta at theta=65, phi=78

The total curve is computed as the sum of the electronic, photonic, bilinear
coupling, and dipole self-energy component differences. Thermal reference
lines are optional and off by default because pQED is used here primarily as a
decomposition model rather than the most quantitative energy benchmark.

Run:
    python plot_pqed_energy_decomposition.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
PQED_DIR = ROOT_DIR / "pqed"

OUTPUT_STEM = "pqed_49_10_energy_decomposition"

HA_TO_KCAL_MOL = 627.509

# Mode-volume conversion used in the related styled QED-CCSD plot.
A0_NM = 0.05291772109
C_NM3 = 4.0 * math.pi * (A0_NM**3)

# Plot toggles.
SHOW_MODE_VOLUME_AXIS = True
SHOW_THERMAL_GUIDES = False
THERMAL_GUIDES_KBT = [2, 4, 6]
RT_KCAL_MOL_298K = 0.00198720425864083 * 298.15

PANEL_SPECS = [
    {
        "file": "pqed_49_10_dir70_31_scan_decomp.csv",
        "isomer_a": "ortho",
        "isomer_b": "meta",
        "theta": 70,
        "phi": 31,
        "title": (
            r"$\Delta E_{\mathrm{ortho-meta}},\; "
            r"\theta=70^{\circ},\; \phi=31^{\circ}$"
        ),
    },
    {
        "file": "pqed_49_10_dir65_78_scan_decomp.csv",
        "isomer_a": "para",
        "isomer_b": "meta",
        "theta": 65,
        "phi": 78,
        "title": (
            r"$\Delta E_{\mathrm{para-meta}},\; "
            r"\theta=65^{\circ},\; \phi=78^{\circ}$"
        ),
    },
]

COMPONENTS = [
    {
        "key": "total",
        "label": "total",
        "color": "gray",
        "linestyle": "-",
        "marker": "o",
        "linewidth": 4.0,
        "markersize": 7.5,
        "zorder": 4,
    },
    {
        "key": "el",
        "label": "electronic",
        "color" :  "black",
        #"color": "#0072B2",
        "linestyle": "-",
        "marker": "s",
        "linewidth": 2.0,
        "markersize": 6.5,
        "zorder": 3,
    },
    {
        "key": "ph",
        "label": "photonic",
        "color": "orange",
        #"color": "#009E73",
        "linestyle": "-",
        "marker": "^",
        "linewidth": 2.0,
        "markersize": 6.5,
        "zorder": 3,
    },
    {
        "key": "blc",
        "label": "bilinear coupling",
        #"color": "#D55E00",
        "color" : "blue",
        "linestyle": "-",
        "marker": "D",
        "linewidth": 2.0,
        "markersize": 6.5,
        "zorder": 3,
    },
    {
        "key": "dse",
        "label": "dipole self energy",
        #"color": "#CC79A7",
        "color" : "red",
        "linestyle": "-",
        "marker": "v",
        "linewidth": 2.2,
        "markersize": 6.5,
        "zorder": 3,
    },
    {
        "key": "blc+dse",
        "label": "blc + dse",
        #"color": "#CC79A7",
        "color" : "purple",
        "linestyle": "--",
        "marker": "v",
        "linewidth": 2.2,
        "markersize": 6.5,
        "zorder": 3,
    },
]


def lambda_to_mode_volume(lambda_au):
    """Convert coupling strength lambda in a.u. to mode volume in nm^3."""
    lambda_au = np.asarray(lambda_au)
    out = np.full_like(lambda_au, np.inf, dtype=float)
    np.divide(C_NM3, lambda_au**2, out=out, where=lambda_au != 0)
    return out


def mode_volume_to_lambda(volume_nm3):
    """Inverse transform for Matplotlib secondary_xaxis."""
    volume_nm3 = np.asarray(volume_nm3)
    ratio = np.full_like(volume_nm3, np.inf, dtype=float)
    np.divide(C_NM3, volume_nm3, out=ratio, where=volume_nm3 != 0)
    return np.sqrt(ratio)


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a decomposition CSV and sort rows by lambda magnitude."""
    if not path.exists():
        raise FileNotFoundError(f"Expected decomposition CSV not found: {path}")

    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    rows.sort(key=lambda row: float(row["lambda_magnitude"]))
    return rows


def component_difference(row: dict[str, str], isomer_a: str, isomer_b: str, component: str) -> float:
    """Return E_A(component) - E_B(component), converted to kcal/mol."""
    value_a = float(row[f"{isomer_a}_E_{component}"])
    value_b = float(row[f"{isomer_b}_E_{component}"])
    return (value_a - value_b) * HA_TO_KCAL_MOL


def build_panel_data(spec: dict) -> dict[str, list[float]]:
    """Compute lambda, component deltas, and total delta for one panel."""
    rows = read_rows(PQED_DIR / spec["file"])
    lambdas = [float(row["lambda_magnitude"]) for row in rows]
    series = {"lambda": lambdas}

    for component in ["el", "ph", "blc", "dse"]:
        series[component] = [
            component_difference(row, spec["isomer_a"], spec["isomer_b"], component)
            for row in rows
        ]

    series["total"] = [
        series["el"][idx] + series["ph"][idx] + series["blc"][idx] + series["dse"][idx]
        for idx in range(len(lambdas))
    ]
    series["blc+dse"] = [
       series["blc"][idx] + series["dse"][idx]
       for idx in range(len(lambdas))
    ]
    return series


def configure_style() -> None:
    """Set plotting defaults for publication-style output."""
    plt.rcParams.update(
        {
            "font.size": 13,
            "font.family": "sans-serif",
            "axes.labelsize": 15,
            "axes.titlesize": 15,
            "axes.linewidth": 1.2,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": True,
            "legend.fontsize": 11,
            "savefig.bbox": "tight",
        }
    )


def add_mode_volume_axis(ax) -> None:
    """Add mode volume as an upper x axis aligned to lambda ticks."""
    secax = ax.secondary_xaxis(
        "top",
        functions=(lambda_to_mode_volume, mode_volume_to_lambda),
    )
    secax.set_xlabel(r"Mode Volume $V$ (nm$^3$)", labelpad=12)

    lambda_ticks = [0.02, 0.04, 0.06, 0.08, 0.10]
    volume_ticks = [lambda_to_mode_volume(value) for value in lambda_ticks]
    secax.set_ticks(volume_ticks)
    secax.set_xticklabels([f"{float(value):.2f}" for value in volume_ticks])


def add_thermal_guides(ax) -> None:
    """Optionally draw thermal stabilization guide lines."""
    for kbt in THERMAL_GUIDES_KBT:
        y = -kbt * RT_KCAL_MOL_298K
        ax.axhline(
            y,
            color="0.25",
            linestyle=":",
            linewidth=1.4,
            alpha=0.55,
            zorder=0,
        )
        ax.text(
            0.021,
            y + 0.25,
            rf"$-{kbt}\;k_B T$",
            fontsize=10,
            color="0.25",
        )


def plot_panel(ax, spec: dict, panel_data: dict[str, list[float]]) -> None:
    """Plot total and component differences for one isomer pair."""
    lambdas = panel_data["lambda"]

    for component in COMPONENTS:
        key = component["key"]
        ax.plot(
            lambdas,
            panel_data[key],
            color=component["color"],
            linestyle=component["linestyle"],
            marker=component["marker"],
            linewidth=component["linewidth"],
            markersize=component["markersize"],
            markerfacecolor=component["color"] if key == "total" else "white",
            markeredgecolor=component["color"],
            markeredgewidth=1.5,
            label=component["label"],
            zorder=component["zorder"],
        )

    ax.axhline(0.0, color="0.45", linestyle="-", linewidth=1.0, alpha=0.5, zorder=0)
    if SHOW_THERMAL_GUIDES:
        add_thermal_guides(ax)
    if SHOW_MODE_VOLUME_AXIS:
        add_mode_volume_axis(ax)

    ax.set_title(spec["title"], pad=14)
    ax.set_xlabel(r"Coupling Strength $\lambda$ (a.u.)", labelpad=8)
    ax.set_xlim(0.015, 0.105)
    ax.set_xticks([0.02, 0.04, 0.06, 0.08, 0.10])
    ax.grid(True, which="both", linestyle=":", alpha=0.45)


def set_shared_y_limits(axes, panel_data_by_spec: list[dict[str, list[float]]]) -> None:
    """Use one y scale across both panels so component magnitudes compare directly."""
    values = []
    for panel_data in panel_data_by_spec:
        for component in COMPONENTS:
            values.extend(panel_data[component["key"]])

    low = min(values)
    high = max(values)
    padding = 0.08 * (high - low)
    y_min = math.floor(low - padding)
    y_max = math.ceil(high + padding)
    for ax in axes:
        ax.set_ylim(y_min, y_max)


def main() -> None:
    configure_style()
    panel_data_by_spec = [build_panel_data(spec) for spec in PANEL_SPECS]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.8), dpi=300, sharey=True)
    set_shared_y_limits(axes, panel_data_by_spec)

    for ax, spec, panel_data in zip(axes, PANEL_SPECS, panel_data_by_spec):
        plot_panel(ax, spec, panel_data)

    axes[0].set_ylabel(r"$\Delta E$ contribution (kcal/mol)", labelpad=10)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(COMPONENTS),
        frameon=True,
        bbox_to_anchor=(0.5, -0.03),
    )

    fig.tight_layout(rect=(0, 0.08, 1, 1))
    png_path = SCRIPT_DIR / f"{OUTPUT_STEM}.png"
    pdf_path = SCRIPT_DIR / f"{OUTPUT_STEM}.pdf"
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Saved {png_path.name} and {pdf_path.name}")


if __name__ == "__main__":
    main()
