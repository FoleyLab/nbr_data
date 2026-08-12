"""
Plot relative isomer energies for the Lambda_Scan_Results dataset.

The script is intentionally data-driven: edit DEFAULT_FIGURES near the bottom
to choose which curves should appear in each figure, then run

    python plot_isomer_energies.py

Figures are written next to this script as both PNG and PDF files.

Supported methods
-----------------
qed_ccsd
    QED-CCSD(2,2) summary CSVs in ../QED_CCSD/summary.
    Required fields: pair, direction, geometry ("relaxed" or "unrelaxed").

qed_dft
    QED-DFT summary CSVs in ../qed_dft_relaxed or ../qed_dft_unrelaxed.
    Required fields: pair, direction, geometry ("relaxed" or "unrelaxed").
    Optional field for relaxed data: zpe=True to use the ZPE-corrected curve.

pqed
    pQED unrelaxed summary CSVs in ../pqed.
    Required fields: pair, direction, Nel, Nph.
    Optional field: CS=True for coherent-state transformed data.

Pair names are "ortho_meta", "para_meta", and, for pQED only, "ortho_para".
Directions are written as tuples, e.g. (70, 31) or (65, 78).
"""

from __future__ import annotations

import csv
from itertools import cycle
from pathlib import Path

import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
OUTPUT_DIR = SCRIPT_DIR

QED_CCSD_SUMMARY_DIR = ROOT_DIR / "QED_CCSD" / "summary"
QED_DFT_RELAXED_DIR = ROOT_DIR / "qed_dft_relaxed"
QED_DFT_UNRELAXED_DIR = ROOT_DIR / "qed_dft_unrelaxed"
PQED_DIR = ROOT_DIR / "pqed"

PAIR_LABELS = {
    "ortho_meta": "ortho - meta",
    "para_meta": "para - meta",
    "ortho_para": "ortho - para",
}

PAIR_FILE_STEMS = {
    "ortho_meta": "ortho-meta",
    "para_meta": "para-meta",
}

PAIR_DFT_UNRELAXED_STEMS = {
    "ortho_meta": "ortho_meta",
    "para_meta": "para_meta",
}

plt.rcParams.update(
    {
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
    }
)

MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


class Series:
    """Simple plottable curve container."""

    def __init__(self, lambda_magnitude: list[float], delta_e: list[float], label: str):
        self.lambda_magnitude = lambda_magnitude
        self.delta_e = delta_e
        self.label = label


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Expected data file not found: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _ordered_series(rows: list[dict[str, str]], value_column: str, label: str) -> Series:
    """Return a sorted, plottable series from a loaded summary table."""
    if not rows:
        raise ValueError("CSV contains no data rows")
    if "lambda_magnitude" not in rows[0]:
        raise KeyError("CSV is missing required column 'lambda_magnitude'")
    if value_column not in rows[0]:
        raise KeyError(f"CSV is missing required column {value_column!r}")

    ordered = sorted(rows, key=lambda row: float(row["lambda_magnitude"]))
    return Series(
        [float(row["lambda_magnitude"]) for row in ordered],
        [float(row[value_column]) for row in ordered],
        label,
    )


def _direction_text(direction: tuple[int, int]) -> str:
    theta, phi = direction
    return f"theta={theta}, phi={phi}"


def _label(entry: dict) -> str:
    if "label" in entry:
        return entry["label"]

    method = entry["method"]
    pair = PAIR_LABELS[entry["pair"]]
    geom = entry.get("geometry")
    direction = _direction_text(entry["direction"])

    if method == "qed_ccsd":
        return f"QED-CCSD(2,2), {geom}, {pair}, {direction}"
    if method == "qed_dft":
        zpe = ", +ZPE" if entry.get("zpe", False) else ""
        return f"QED-DFT, {geom}{zpe}, {pair}, {direction}"
    if method == "pqed":
        cs = ", CS" if entry.get("CS", False) else ""
        return f"pQED ({entry['Nel']}e,{entry['Nph']}ph{cs}), {pair}, {direction}"
    raise ValueError(f"Unknown method {method!r}")


def load_qed_ccsd(entry: dict) -> Series:
    pair = entry["pair"]
    theta, phi = entry["direction"]
    geometry = entry["geometry"]
    path = QED_CCSD_SUMMARY_DIR / f"{geometry}_dir_{theta}_{phi}_QEDCCSD22.csv"
    column = f"dE_{pair}_kcal/mol"
    return _ordered_series(_read_csv(path), column, _label(entry))


def load_qed_dft(entry: dict) -> Series:
    pair = entry["pair"]
    theta, phi = entry["direction"]
    geometry = entry["geometry"]

    if pair not in PAIR_FILE_STEMS:
        raise ValueError("QED-DFT data only contains ortho_meta and para_meta pairs")

    if geometry == "relaxed":
        path = QED_DFT_RELAXED_DIR / f"{PAIR_FILE_STEMS[pair]}_{theta}_{phi}_hartree.csv"
        suffix = "zpe" if entry.get("zpe", False) else "raw"
        column = f"dE_{pair}_{suffix}_kcal/mol"
    elif geometry == "unrelaxed":
        stem = PAIR_DFT_UNRELAXED_STEMS[pair]
        path = QED_DFT_UNRELAXED_DIR / f"{stem}_dir{theta}_{phi}_scan_qed_dft_no_relax.csv"
        column = f"dE_{pair}_kcal/mol"
    else:
        raise ValueError("QED-DFT geometry must be 'relaxed' or 'unrelaxed'")

    return _ordered_series(_read_csv(path), column, _label(entry))


def load_pqed(entry: dict) -> Series:
    pair = entry["pair"]
    theta, phi = entry["direction"]
    cs_suffix = "_CS" if entry.get("CS", False) else ""
    path = PQED_DIR / f"pqed_{entry['Nel']}_{entry['Nph']}_dir{theta}_{phi}_scan{cs_suffix}.csv"
    column = f"dE_{pair}_kcal/mol"
    return _ordered_series(_read_csv(path), column, _label(entry))


LOADERS = {
    "qed_ccsd": load_qed_ccsd,
    "qed_dft": load_qed_dft,
    "pqed": load_pqed,
}


def load_series(entry: dict) -> Series:
    method = entry["method"]
    if method not in LOADERS:
        raise ValueError(f"Unknown method {method!r}; choose from {sorted(LOADERS)}")
    return LOADERS[method](entry)


def make_plot(output_name: str, series_entries: list[dict], title: str | None = None) -> None:
    """Overlay configured series and save one PNG plus one PDF."""
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    marker_cycle = cycle(MARKERS)

    for entry in series_entries:
        series = load_series(entry)
        ax.plot(
            series.lambda_magnitude,
            series.delta_e,
            marker=next(marker_cycle),
            label=series.label,
        )

    ax.axhline(0.0, color="0.55", linewidth=0.9, linestyle="--", zorder=0)
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


# Edit this mapping to add or remove figures. Each key is an output filename
# stem; each value is a list of curve definitions.
DEFAULT_FIGURES = {
    "01_relaxed_vs_unrelaxed_no_zpe": [
        dict(method="qed_dft", geometry="relaxed", pair="ortho_meta", direction=(70, 31)),
        dict(method="qed_dft", geometry="unrelaxed", pair="ortho_meta", direction=(70, 31)),
        dict(method="qed_dft", geometry="relaxed", pair="para_meta", direction=(65, 78)),
        dict(method="qed_dft", geometry="unrelaxed", pair="para_meta", direction=(65, 78)),
    ],
    "02_relaxed_zpe_vs_no_zpe": [
        dict(method="qed_dft", geometry="relaxed", pair="ortho_meta", direction=(70, 31), zpe=False),
        dict(method="qed_dft", geometry="relaxed", pair="ortho_meta", direction=(70, 31), zpe=True),
        dict(method="qed_dft", geometry="relaxed", pair="para_meta", direction=(65, 78), zpe=False),
        dict(method="qed_dft", geometry="relaxed", pair="para_meta", direction=(65, 78), zpe=True),
    ],
    "03_qeddft_vs_pqed_ortho_meta": [
        dict(method="qed_dft", geometry="relaxed", pair="ortho_meta", direction=(70, 31), zpe=True),
        dict(method="qed_dft", geometry="unrelaxed", pair="ortho_meta", direction=(70, 31)),
        dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=10),
        dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=3),
    ],
    "04_qeddft_vs_pqed_para_meta": [
        dict(method="qed_dft", geometry="relaxed", pair="para_meta", direction=(65, 78), zpe=True),
        dict(method="qed_dft", geometry="unrelaxed", pair="para_meta", direction=(65, 78)),
        dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=10),
        dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=3),
    ],
    "05_qed_ccsd_relaxed_vs_unrelaxed": [
        dict(method="qed_ccsd", geometry="relaxed", pair="ortho_meta", direction=(70, 31)),
        dict(method="qed_ccsd", geometry="unrelaxed", pair="ortho_meta", direction=(70, 31)),
        dict(method="qed_ccsd", geometry="relaxed", pair="para_meta", direction=(65, 78)),
        dict(method="qed_ccsd", geometry="unrelaxed", pair="para_meta", direction=(65, 78)),
    ],
    "05a_relaxed_vs_unrelaxed_ortho_meta": [
        dict(method="qed_dft", geometry="unrelaxed", pair="ortho_meta", direction=(70, 31), zpe=False),
        dict(method="qed_dft", geometry="relaxed", pair="ortho_meta", direction=(70, 31), zpe=False),
        dict(method="qed_ccsd", geometry="relaxed", pair="ortho_meta", direction=(70, 31)),
        dict(method="qed_ccsd", geometry="unrelaxed", pair="ortho_meta", direction=(70, 31)),
    ],
    "05b_relaxed_vs_unrelaxed_ortho_meta": [
        dict(method="qed_dft", geometry="unrelaxed", pair="para_meta", direction=(65, 78), zpe=False),
        dict(method="qed_dft", geometry="relaxed", pair="para_meta", direction=(65, 78), zpe=False),
        dict(method="qed_ccsd", geometry="relaxed", pair="para_meta", direction=(65, 78)),
        dict(method="qed_ccsd", geometry="unrelaxed", pair="para_meta", direction=(65, 78)),
    ],

    "06_pqed_vs_qed_ccsd_ortho_meta": [
        dict(method="qed_ccsd", geometry="relaxed", pair="ortho_meta", direction=(70, 31)),
        dict(method="qed_ccsd", geometry="unrelaxed", pair="ortho_meta", direction=(70, 31)),
        dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=10),
        dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=3),
    ],
    "07_pqed_vs_qed_ccsd_para_meta": [
        dict(method="qed_ccsd", geometry="relaxed", pair="para_meta", direction=(65, 78)),
        dict(method="qed_ccsd", geometry="unrelaxed", pair="para_meta", direction=(65, 78)),
        dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=10),
        dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=3),
    ],
    "08_pqed_coherent_state_check": [
        dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=10, CS=False),
        dict(method="pqed", pair="ortho_meta", direction=(70, 31), Nel=49, Nph=10, CS=True),
        dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=10, CS=False),
        dict(method="pqed", pair="para_meta", direction=(65, 78), Nel=49, Nph=10, CS=True),
    ],
}


def main() -> None:
    for output_name, entries in DEFAULT_FIGURES.items():
        make_plot(output_name, entries)


if __name__ == "__main__":
    main()
