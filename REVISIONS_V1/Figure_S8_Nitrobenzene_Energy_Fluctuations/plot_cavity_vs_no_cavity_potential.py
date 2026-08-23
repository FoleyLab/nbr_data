"""Compare cavity and cavity-free potential energies along one QED-DFT trajectory.

The CSV is assumed to contain cavity and cavity-free single-point energies evaluated
at the same sampled nuclear geometries R_t.  The resulting figure is designed for
Supporting Information and directly addresses the time-averaged energetic shift
induced by the cavity.

Expected columns:
    frame
    cavity_energy_hartree
    no_cavity_energy_hartree
    delta_cavity_minus_no_cavity_hartree
    delta_cavity_minus_no_cavity_kcal_mol

Outputs:
    cavity_vs_no_cavity_potential_energy.png
    cavity_vs_no_cavity_potential_energy.pdf
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_CSV = SCRIPT_DIR / "nitrobenzene_direction_A_wb97x_d_cavity_free_energies.csv"
OUTPUT_STEM = "cavity_vs_no_cavity_potential_energy"

HARTREE_TO_KCAL_MOL = 627.5094740631

# If the MD time step is known, set this to the time in fs between *original*
# trajectory frames, e.g. 0.5.  Leave as None to plot trajectory frame number.
DT_FS_PER_FRAME: float | None = None


def read_data(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    return {
        "frame": np.asarray([float(row["frame"]) for row in rows]),
        "cavity": np.asarray([float(row["cavity_energy_hartree"]) for row in rows]),
        "nocavity": np.asarray([float(row["no_cavity_energy_hartree"]) for row in rows]),
        "delta_kcal": np.asarray(
            [float(row["delta_cavity_minus_no_cavity_kcal_mol"]) for row in rows]
        ),
    }


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 12,
            "font.family": "sans-serif",
            "axes.labelsize": 13,
            "axes.linewidth": 1.1,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "ytick.right": True,
            "legend.fontsize": 10,
            "savefig.bbox": "tight",
        }
    )


def main() -> None:
    data = read_data(INPUT_CSV)
    frame = data["frame"]
    cavity = data["cavity"]
    nocavity = data["nocavity"]
    delta = data["delta_kcal"]

    if DT_FS_PER_FRAME is None:
        x = frame
        xlabel = "Trajectory frame"
    else:
        x = frame * DT_FS_PER_FRAME
        xlabel = "Time (fs)"

    # Use one common reference for both curves so their vertical separation remains
    # physically meaningful.  Choosing the time-averaged cavity-free energy makes
    # the cavity-free mean exactly zero by construction.
    e_ref = np.mean(nocavity)
    cavity_rel = (cavity - e_ref) * HARTREE_TO_KCAL_MOL
    nocavity_rel = (nocavity - e_ref) * HARTREE_TO_KCAL_MOL

    mean_cavity_rel = float(np.mean(cavity_rel))
    mean_nocavity_rel = float(np.mean(nocavity_rel))
    mean_delta = float(np.mean(delta))
    std_delta = float(np.std(delta, ddof=1))

    configure_style()
    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(9.0, 7.0),
        dpi=300,
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [1.15, 1.0]},
    )

    # Panel a: direct comparison at identical R_t.
    line_nocav, = ax_top.plot(x, nocavity_rel, linewidth=1.4, label="Cavity-free")
    line_cav, = ax_top.plot(x, cavity_rel, linewidth=1.4, label="Cavity")
    ax_top.axhline(
        mean_nocavity_rel,
        linestyle="--",
        linewidth=1.2,
        color=line_nocav.get_color(),
    )
    ax_top.axhline(
        mean_cavity_rel,
        linestyle="--",
        linewidth=1.2,
        color=line_cav.get_color(),
    )
    ax_top.set_ylabel(
        r"$E(\mathbf{R}_t)-\langle E_{\mathrm{no\ cav}}\rangle$"
        "\n(kcal mol$^{-1}$)"
    )
    ax_top.legend(loc="best", frameon=False)
    ax_top.text(
        0.02,
        0.22,
        "Same nuclear geometries in both calculations",
        transform=ax_top.transAxes,
        va="top",
        fontsize=10,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8},
    )

    # Panel b: cavity-induced shift itself.
    line_delta, = ax_bottom.plot(x, delta, linewidth=1.5)
    ax_bottom.axhline(
        mean_delta,
        linestyle="--",
        linewidth=1.4,
        color=line_delta.get_color(),
        label=rf"Mean = {mean_delta:.1f} kcal mol$^{{-1}}$",
    )
    ax_bottom.fill_between(
        x,
        mean_delta - std_delta,
        mean_delta + std_delta,
        color=line_delta.get_color(),
        alpha=0.12,
        linewidth=0,
        label=rf"$\pm 1\sigma$ = {std_delta:.1f} kcal mol$^{{-1}}$",
    )
    ax_bottom.set_xlabel(xlabel)
    ax_bottom.set_ylabel(
        r"$E_{\mathrm{cav}}(\mathbf{R}_t)-E_{\mathrm{no\ cav}}(\mathbf{R}_t)$"
        "\n(kcal mol$^{-1}$)"
    )
    ax_bottom.legend(loc="best", frameon=False)

    for label, ax in zip(("(a)", "(b)"), (ax_top, ax_bottom)):
        ax.text(0.01, 0.98, label, transform=ax.transAxes, va="top", fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.35)


    png_path = SCRIPT_DIR / f"{OUTPUT_STEM}.png"
    pdf_path = SCRIPT_DIR / f"{OUTPUT_STEM}.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)

    print(f"N sampled geometries: {len(frame)}")
    print(f"Mean cavity-free energy: {np.mean(nocavity):.12f} Eh")
    print(f"Mean cavity energy:      {np.mean(cavity):.12f} Eh")
    print(f"Mean cavity shift:       {mean_delta:.3f} kcal/mol")
    print(f"Std. dev. cavity shift:  {std_delta:.3f} kcal/mol")
    print(f"Range of cavity shift:   {np.min(delta):.3f} to {np.max(delta):.3f} kcal/mol")
    print(f"Saved {png_path.name} and {pdf_path.name}")


if __name__ == "__main__":
    main()
