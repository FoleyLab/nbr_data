#!/usr/bin/env python
"""
analyze.py

Analyze the cavity opt+frequency campaign: relative energies of two isomer
pairs versus cavity coupling magnitude, along one cavity orientation each.

Curves produced (relative energy = A - B):
  (70, 31)  ortho - meta   raw           -> Line 1
  (70, 31)  ortho - meta   ZPE-corrected -> Line 2
  (65, 78)  para  - meta   raw           -> Line 3
  (65, 78)  para  - meta   ZPE-corrected -> Line 4

"raw" uses opt_status.json:final_energy_hartree.
"ZPE-corrected" adds frequencies.json:zpe_hartree to each isomer's energy
before differencing: (E_A + zpe_A) - (E_B + zpe_B).

Outputs (written to analysis/qed_dft_relaxed/):
  * <pair>_<orientation>_hartree.csv  -- full per-lambda table in Hartree
  * relative_energies_kcal.png        -- the four curves in kcal/mol
  * availability.md                   -- which lambda points / cells are present

The plot is in kcal/mol (Hartree * 627.509); the tables are in Hartree, as
requested. Missing cells are skipped gracefully and reported, never invented.
"""

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


HARTREE_TO_KCAL = 627.509

BASE = Path(__file__).resolve().parent
RUNS = BASE / "runs"
OUT = BASE / "analysis" / "qed_dft_relaxed"

MAGNITUDES = [0.02, 0.04, 0.06, 0.08, 0.10]

# (label, isomer_A, isomer_B, theta, phi)
PAIRS = [
    ("ortho-meta", "ortho", "meta", 70, 31),
    ("para-meta", "para", "meta", 65, 78),
]


def cell_id(isomer, theta, phi, mag):
    return f"{isomer}_{theta}_{phi}_lam{mag:.2f}"


def load_cell(isomer, theta, phi, mag):
    """Return dict with energy (Ha), zpe (Ha), converged, n_imag; None where absent."""
    d = RUNS / cell_id(isomer, theta, phi, mag)
    out = {"E": None, "zpe": None, "converged": None, "n_imag": None,
           "has_opt": False, "has_freq": False}

    opt_path = d / "opt_status.json"
    if opt_path.exists():
        opt = json.loads(opt_path.read_text())
        out["E"] = opt.get("final_energy_hartree")
        out["converged"] = opt.get("converged")
        out["has_opt"] = True

    freq_path = d / "frequencies.json"
    if freq_path.exists():
        freq = json.loads(freq_path.read_text())
        out["zpe"] = freq.get("zpe_hartree")
        out["n_imag"] = freq.get("n_imaginary")
        out["has_freq"] = True

    return out


def build_pair_table(isomer_A, isomer_B, theta, phi):
    """Assemble the per-lambda table (Hartree) for one isomer pair."""
    rows = []
    for mag in MAGNITUDES:
        a = load_cell(isomer_A, theta, phi, mag)
        b = load_cell(isomer_B, theta, phi, mag)

        dE_raw = None
        if a["E"] is not None and b["E"] is not None:
            dE_raw = a["E"] - b["E"]

        dE_zpe = None
        if None not in (a["E"], b["E"], a["zpe"], b["zpe"]):
            dE_zpe = (a["E"] + a["zpe"]) - (b["E"] + b["zpe"])

        rows.append({
            "magnitude": mag,
            "E_A_Ha": a["E"], "E_B_Ha": b["E"],
            "zpe_A_Ha": a["zpe"], "zpe_B_Ha": b["zpe"],
            "dE_raw_Ha": dE_raw, "dE_zpe_Ha": dE_zpe,
            "A_converged": a["converged"], "B_converged": b["converged"],
            "A_has_freq": a["has_freq"], "B_has_freq": b["has_freq"],
        })
    return rows


def write_csv(path, label, isomer_A, isomer_B, theta, phi, rows):
    cols = ["magnitude", "E_A_Ha", "E_B_Ha", "zpe_A_Ha", "zpe_B_Ha",
            "dE_raw_Ha", "dE_zpe_Ha", "A_converged", "B_converged",
            "A_has_freq", "B_has_freq"]
    header_note = f"# {label}  ({isomer_A} - {isomer_B}) at (theta={theta}, phi={phi}); energies in Hartree"
    with open(path, "w") as f:
        f.write(header_note + "\n")
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")


def make_plot(all_tables, path):
    plt.figure(figsize=(8, 6))

    styles = {
        "ortho-meta": {"color": "#1f77b4"},
        "para-meta": {"color": "#d62728"},
    }

    for (label, iA, iB, theta, phi), rows in all_tables:
        color = styles[label]["color"]

        def both_converged(r):
            return bool(r["A_converged"]) and bool(r["B_converged"])

        # ---- raw curve ----
        raw_pts = [(r["magnitude"], r["dE_raw_Ha"] * HARTREE_TO_KCAL, both_converged(r))
                   for r in rows if r["dE_raw_Ha"] is not None]
        if raw_pts:
            xs = [p[0] for p in raw_pts]
            ys = [p[1] for p in raw_pts]
            plt.plot(xs, ys, "-", color=color, label=f"{label} raw ({theta},{phi})")
            # filled markers = both optimizations converged; hollow = at least one not.
            solid = [(x, y) for x, y, ok in raw_pts if ok]
            shaky = [(x, y) for x, y, ok in raw_pts if not ok]
            if solid:
                plt.plot(*zip(*solid), "o", color=color)
            if shaky:
                plt.plot(*zip(*shaky), "o", markerfacecolor="none",
                         markeredgecolor=color)

        # ---- ZPE-corrected curve ----
        zpe_pts = [(r["magnitude"], r["dE_zpe_Ha"] * HARTREE_TO_KCAL, both_converged(r))
                   for r in rows if r["dE_zpe_Ha"] is not None]
        if zpe_pts:
            xs = [p[0] for p in zpe_pts]
            ys = [p[1] for p in zpe_pts]
            plt.plot(xs, ys, "--", color=color, alpha=0.7,
                     label=f"{label} +ZPE ({theta},{phi})")
            solid = [(x, y) for x, y, ok in zpe_pts if ok]
            shaky = [(x, y) for x, y, ok in zpe_pts if not ok]
            if solid:
                plt.plot(*zip(*solid), "s", color=color, alpha=0.7)
            if shaky:
                plt.plot(*zip(*shaky), "s", markerfacecolor="none",
                         markeredgecolor=color, alpha=0.7)

    plt.xlabel(r"Cavity coupling magnitude $|\lambda|$ (a.u.)")
    plt.ylabel(r"Relative energy $A-B$ (kcal/mol)")
    plt.title("Isomer relative energies vs cavity coupling")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.figtext(0.5, 0.005,
                "Filled markers: both optimizations converged (projected |g| < 5e-4).  "
                "Hollow markers: >=1 optimization NOT converged.",
                ha="center", fontsize=8, style="italic")
    plt.tight_layout(rect=(0, 0.03, 1, 1))
    plt.savefig(path, dpi=150)
    plt.close()


def write_availability(all_tables, path):
    lines = ["# Data availability", ""]
    for (label, iA, iB, theta, phi), rows in all_tables:
        lines.append(f"## {label}  ({iA} - {iB}) at (theta={theta}, phi={phi})")
        lines.append("")
        lines.append("| |lambda| | raw (A-B) | +ZPE (A-B) | missing for this point |")
        lines.append("| --- | --- | --- | --- |")
        for r in rows:
            missing = []
            if r["E_A_Ha"] is None:
                missing.append(f"{iA} opt")
            if r["E_B_Ha"] is None:
                missing.append(f"{iB} opt")
            if r["zpe_A_Ha"] is None:
                missing.append(f"{iA} freq")
            if r["zpe_B_Ha"] is None:
                missing.append(f"{iB} freq")
            raw_ok = "yes" if r["dE_raw_Ha"] is not None else "NO"
            zpe_ok = "yes" if r["dE_zpe_Ha"] is not None else "NO"
            lines.append(f"| {r['magnitude']:.2f} | {raw_ok} | {zpe_ok} | "
                         f"{', '.join(missing) if missing else '-'} |")
        lines.append("")
    path.write_text("\n".join(lines))


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    all_tables = []
    for label, iA, iB, theta, phi in PAIRS:
        rows = build_pair_table(iA, iB, theta, phi)
        all_tables.append(((label, iA, iB, theta, phi), rows))
        csv_path = OUT / f"{label}_{theta}_{phi}_hartree.csv"
        write_csv(csv_path, label, iA, iB, theta, phi, rows)

    make_plot(all_tables, OUT / "relative_energies_kcal.png")
    write_availability(all_tables, OUT / "availability.md")

    # Console summary
    print("=== Curve availability (points usable / 5) ===")
    for (label, iA, iB, theta, phi), rows in all_tables:
        raw_n = sum(1 for r in rows if r["dE_raw_Ha"] is not None)
        zpe_n = sum(1 for r in rows if r["dE_zpe_Ha"] is not None)
        print(f"  {label} ({theta},{phi}):  raw {raw_n}/5   +ZPE {zpe_n}/5")
    print(f"\nWrote outputs to: {OUT}")


if __name__ == "__main__":
    main()
