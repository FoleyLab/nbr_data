#!/usr/bin/env python3
"""
Direct-computation pQED energy-decomposition plot.

This script makes the same style of component-decomposition plot as
plot_pqed_decomp.py, but computes the pQED component expectation values
directly from the EOM-CCSD data instead of reading precomputed decomposition
CSVs.  It can plot the non-CS Hamiltonian, the coherent-state transformed
Hamiltonian, or both as separate output files.

Run:
    conda run -n p4dev python plot_pqed_decomp_direct.py
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "matplotlib-pqed-direct"),
)

import matplotlib.pyplot as plt
import numpy as np


sys.path.insert(0, str(SCRIPT_DIR))

import pf_isomer_scan as pfscan  # noqa: E402


# ============================================================================
#                              CONFIGURATION
# ============================================================================

# --- Isomer pair and cavity-field direction ---------------------------------
ISOMER_A = "para"
ISOMER_B = "meta"
THETA = 65.0
PHI = 78.0

# --- PF Hamiltonian size knobs ----------------------------------------------
NUM_ELECTRONIC_STATES = 49
NUM_FOCK_STATES = 2

# --- Photon frequency (Hartree), fit from the EOM-CCSD data -----------------
OMEGA = 0.066148

# --- Lambda magnitudes to scan over -----------------------------------------
LAMBDA_MAGNITUDES = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]

# --- Plot mode: "non-cs", "cs", or "both" ---------------------------------
TRANSFORMATION = "both"

# --- Shift total/electronic curves to their lambda=0 baseline ---------------
SUBTRACT_ZERO_BASELINE = True

# --- Output naming -----------------------------------------------------------
OUTPUT_STEM = "pqed_direct_decomposition"
SAVE_PNG = True
SAVE_PDF = True

# --- Plot toggles ------------------------------------------------------------
SHOW_MODE_VOLUME_AXIS = True
SHOW_THERMAL_GUIDES = False
THERMAL_GUIDES_KBT = [2, 4, 6]


# ============================================================================
#                              CONSTANTS / STYLE
# ============================================================================

HA_TO_KCAL_MOL = 627.509
RT_KCAL_MOL_298K = 0.00198720425864083 * 298.15

# Mode-volume conversion used in the related styled QED-CCSD plot.
A0_NM = 0.05291772109
C_NM3 = 4.0 * math.pi * (A0_NM**3)

ISOMER_LABELS = {
    "ortho": "Ortho",
    "meta": "Meta",
    "para": "Para",
}

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
        "color": "black",
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
        "linestyle": "-",
        "marker": "^",
        "linewidth": 2.0,
        "markersize": 6.5,
        "zorder": 3,
    },
    {
        "key": "blc",
        "label": "bilinear coupling",
        "color": "blue",
        "linestyle": "-",
        "marker": "D",
        "linewidth": 2.0,
        "markersize": 6.5,
        "zorder": 3,
    },
    {
        "key": "dse",
        "label": "dipole self energy",
        "color": "red",
        "linestyle": "-",
        "marker": "v",
        "linewidth": 2.2,
        "markersize": 6.5,
        "zorder": 3,
    },
    {
        "key": "blc+dse",
        "label": "blc + dse",
        "color": "purple",
        "linestyle": "--",
        "marker": "v",
        "linewidth": 2.2,
        "markersize": 6.5,
        "zorder": 3,
    },
]


# ============================================================================
#                              DATA / COMPUTATION
# ============================================================================

def normalize_isomer_name(name: str) -> str:
    """Convert user-facing isomer names to the labels used by pf_isomer_scan."""
    key = name.strip().lower()
    try:
        return ISOMER_LABELS[key]
    except KeyError as exc:
        valid = ", ".join(sorted(ISOMER_LABELS))
        raise ValueError(f"Unknown isomer {name!r}. Expected one of: {valid}") from exc


def load_isomer_data(num_electronic_states: int) -> dict[str, dict[str, np.ndarray]]:
    """Parse and symmetrize EOM-CCSD data for each isomer."""
    isomer_data = {}
    for label, fname in pfscan.ISOMERS.items():
        path = SCRIPT_DIR / pfscan.DATA_DIR / fname
        ref_E, corr_E, excit_E, electronic_dipole, total_dipoles = pfscan.parse_cq_h5_data(
            str(path),
            verbose=True,
        )

        electronic_dipoles_sym = 0.5 * (
            electronic_dipole + np.transpose(electronic_dipole, axes=(1, 0, 2))
        )
        total_dipoles_sym = 0.5 * (
            total_dipoles + np.transpose(total_dipoles, axes=(1, 0, 2))
        )

        e_el = pfscan.build_electronic_energies(ref_E, corr_E, excit_E)
        if num_electronic_states > len(e_el):
            raise ValueError(
                f"{label}: NUM_ELECTRONIC_STATES={num_electronic_states} "
                f"exceeds available states ({len(e_el)})."
            )

        isomer_data[label] = {
            "E_el": e_el,
            "electronic_dipoles": electronic_dipoles_sym,
            "total_dipoles": total_dipoles_sym,
        }

    return isomer_data


def compute_components_for_isomer(
    data: dict[str, np.ndarray],
    lambda_vec: np.ndarray,
    transformation: str,
) -> dict[str, float]:
    """Diagonalize the selected PF Hamiltonian and return component energies."""
    (
        h_total,
        h_total_cs,
        h_el_full,
        h_ph_full,
        h_blc,
        h_dse,
        h_blc_cs,
        h_dse_cs,
    ) = pfscan.build_PF_Hamiltonian(
        dim_ph=NUM_FOCK_STATES,
        dim_el=NUM_ELECTRONIC_STATES,
        omega=OMEGA,
        lambda_vec=lambda_vec,
        e_el=data["E_el"],
        electronic_dipoles=data["electronic_dipoles"],
        total_dipoles=data["total_dipoles"],
    )

    if transformation == "cs":
        _, coeffs = np.linalg.eigh(h_total_cs)
        psi0 = coeffs[:, 0]
        h_blc_selected = h_blc_cs
        h_dse_selected = h_dse_cs
    elif transformation == "non-cs":
        _, coeffs = np.linalg.eigh(h_total)
        psi0 = coeffs[:, 0]
        h_blc_selected = h_blc
        h_dse_selected = h_dse
    else:
        raise ValueError(f"Unknown transformation: {transformation!r}")

    return {
        "el": float(np.real(psi0.conj().T @ h_el_full @ psi0)),
        "ph": float(np.real(psi0.conj().T @ h_ph_full @ psi0)),
        "blc": float(np.real(psi0.conj().T @ h_blc_selected @ psi0)),
        "dse": float(np.real(psi0.conj().T @ h_dse_selected @ psi0)),
    }


def build_panel_data(
    isomer_data: dict[str, dict[str, np.ndarray]],
    isomer_a: str,
    isomer_b: str,
    transformation: str,
) -> dict[str, list[float]]:
    """Compute lambda, component deltas, and total delta for one plot."""
    lambda_dir = pfscan.generate_lambda_vec_from_theta_and_phi(THETA, PHI)
    lambdas = sorted(float(value) for value in LAMBDA_MAGNITUDES)
    series = {
        "lambda": lambdas,
        "el": [],
        "ph": [],
        "blc": [],
        "dse": [],
        "blc+dse": [],
        "total": [],
    }

    zero_baseline = {"el": 0.0, "total": 0.0}
    if SUBTRACT_ZERO_BASELINE:
        zero_vec = lambda_dir * 0.0
        comp_a_zero = compute_components_for_isomer(
            isomer_data[isomer_a],
            zero_vec,
            transformation,
        )
        comp_b_zero = compute_components_for_isomer(
            isomer_data[isomer_b],
            zero_vec,
            transformation,
        )
        zero_baseline["el"] = (
            comp_a_zero["el"] - comp_b_zero["el"]
        ) * HA_TO_KCAL_MOL
        zero_baseline["total"] = sum(
            (comp_a_zero[key] - comp_b_zero[key]) * HA_TO_KCAL_MOL
            for key in ["el", "ph", "blc", "dse"]
        )

    for lambda_magnitude in lambdas:
        lambda_vec = lambda_dir * lambda_magnitude
        comp_a = compute_components_for_isomer(
            isomer_data[isomer_a],
            lambda_vec,
            transformation,
        )
        comp_b = compute_components_for_isomer(
            isomer_data[isomer_b],
            lambda_vec,
            transformation,
        )

        for key in ["el", "ph", "blc", "dse"]:
            series[key].append((comp_a[key] - comp_b[key]) * HA_TO_KCAL_MOL)

        idx = len(series["total"])
        series["total"].append(
            series["el"][idx]
            + series["ph"][idx]
            + series["blc"][idx]
            + series["dse"][idx]
        )
        series["blc+dse"].append(series["blc"][idx] + series["dse"][idx])

        if SUBTRACT_ZERO_BASELINE:
            series["el"][idx] -= zero_baseline["el"]
            series["total"][idx] -= zero_baseline["total"]

        print(
            f"  {transformation:6s} lambda={lambda_magnitude:.4f}: "
            f"total dE = {series['total'][idx]: .8f} kcal/mol"
        )

    return series


# ============================================================================
#                              PLOTTING
# ============================================================================

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

    volume_ticks = [lambda_to_mode_volume(value) for value in LAMBDA_MAGNITUDES]
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
            min(LAMBDA_MAGNITUDES) * 1.05,
            y + 0.25,
            rf"$-{kbt}\;k_B T$",
            fontsize=10,
            color="0.25",
        )


def plot_panel(ax, title: str, panel_data: dict[str, list[float]]) -> None:
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

    x_padding = 0.05 * (max(lambdas) - min(lambdas))
    ax.set_title(title, pad=14)
    ax.set_xlabel(r"Coupling Strength $\lambda$ (a.u.)", labelpad=8)
    ax.set_ylabel(r"$\Delta E$ contribution (kcal/mol)", labelpad=10)
    ax.set_xlim(min(lambdas) - x_padding, max(lambdas) + x_padding)
    ax.set_xticks(lambdas)
    ax.grid(True, which="both", linestyle=":", alpha=0.45)


def title_for_plot(isomer_a: str, isomer_b: str, transformation: str) -> str:
    """Build the mathtext title used above the plot."""
    label_a = isomer_a.lower()
    label_b = isomer_b.lower()
    suffix = "CS" if transformation == "cs" else "PF"
    return (
        rf"$\Delta E_{{\mathrm{{{label_a}-{label_b}}}}},\; "
        rf"\theta={THETA:.0f}^{{\circ}},\; \phi={PHI:.0f}^{{\circ}}$"
        f" ({suffix})"
    )


def output_paths(isomer_a: str, isomer_b: str, transformation: str) -> tuple[Path, Path]:
    """Return PNG and PDF output paths for one plot."""
    theta_tag = f"{THETA:g}".replace(".", "p")
    phi_tag = f"{PHI:g}".replace(".", "p")
    suffix = (
        f"{isomer_a.lower()}_{isomer_b.lower()}"
        f"_dir{theta_tag}_{phi_tag}"
        f"_{NUM_ELECTRONIC_STATES}_{NUM_FOCK_STATES}"
        f"_{transformation.replace('-', '')}"
    )
    stem = f"{OUTPUT_STEM}_{suffix}"
    return SCRIPT_DIR / f"{stem}.png", SCRIPT_DIR / f"{stem}.pdf"


def save_plot(
    panel_data: dict[str, list[float]],
    isomer_a: str,
    isomer_b: str,
    transformation: str,
) -> None:
    """Create and save one direct-computation decomposition plot."""
    configure_style()
    fig, ax = plt.subplots(1, 1, figsize=(7.4, 6.8), dpi=300)
    plot_panel(ax, title_for_plot(isomer_a, isomer_b, transformation), panel_data)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        frameon=True,
        bbox_to_anchor=(0.5, -0.03),
    )

    fig.tight_layout(rect=(0, 0.10, 1, 1))
    png_path, pdf_path = output_paths(isomer_a, isomer_b, transformation)
    if SAVE_PNG:
        fig.savefig(png_path)
    if SAVE_PDF:
        fig.savefig(pdf_path)
    plt.close(fig)

    saved = []
    if SAVE_PNG:
        saved.append(png_path.name)
    if SAVE_PDF:
        saved.append(pdf_path.name)
    if saved:
        print(f"Saved {', '.join(saved)}")


# ============================================================================
#                                 MAIN
# ============================================================================

def transformations_to_run() -> list[str]:
    """Normalize TRANSFORMATION into the list of plot modes to generate."""
    key = TRANSFORMATION.strip().lower()
    if key == "both":
        return ["non-cs", "cs"]
    if key in {"non-cs", "cs"}:
        return [key]
    raise ValueError('TRANSFORMATION must be "non-cs", "cs", or "both"')


def main() -> None:
    isomer_a = normalize_isomer_name(ISOMER_A)
    isomer_b = normalize_isomer_name(ISOMER_B)
    if isomer_a == isomer_b:
        raise ValueError("ISOMER_A and ISOMER_B must be different.")

    print(
        f"Computing direct pQED decomposition for {isomer_a.lower()} - {isomer_b.lower()} "
        f"at theta={THETA:g}, phi={PHI:g}"
    )
    print(
        f"  electronic states = {NUM_ELECTRONIC_STATES}, "
        f"Fock states = {NUM_FOCK_STATES}, omega = {OMEGA}"
    )
    print(f"  lambda magnitudes = {LAMBDA_MAGNITUDES}")

    isomer_data = load_isomer_data(NUM_ELECTRONIC_STATES)
    for transformation in transformations_to_run():
        panel_data = build_panel_data(
            isomer_data=isomer_data,
            isomer_a=isomer_a,
            isomer_b=isomer_b,
            transformation=transformation,
        )
        save_plot(panel_data, isomer_a, isomer_b, transformation)


if __name__ == "__main__":
    main()
