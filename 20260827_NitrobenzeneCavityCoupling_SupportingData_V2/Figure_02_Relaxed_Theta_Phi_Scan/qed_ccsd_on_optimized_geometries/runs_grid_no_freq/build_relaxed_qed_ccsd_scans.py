#!/usr/bin/env python3
"""Parse SCF/CCSD energies from optimized-geometry QED-CCSD runs and build
relaxed_qed_ccsd_intermediate_scans.csv to match unrelaxed_..._scans.csv."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

ISOMERS = ("ortho", "meta", "para")

SCF_RE = re.compile(r"\*\* Total SCF energy\s*=\s*([-0-9.]+)")
CCSD_RE = re.compile(r"CCSD total energy\s*/\s*hartree\s*=\s*([-0-9.]+)")

FIELD_NAMES = [
    "theta",
    "phi",
    "E_SCF_ortho_int",
    "E_CCSD_ortho_int",
    "E_SCF_meta_int",
    "E_CCSD_meta_int",
    "E_SCF_para_int",
    "E_CCSD_para_int",
]


def parse_runs_dir(runs_dir: Path, grid_path: Path) -> list[dict[str, object]]:
    with grid_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        grid = [(float(row["theta"]), float(row["phi"])) for row in reader]

    rows = []
    for theta, phi in grid:
        row: dict[str, object] = {"theta": f"{theta:.1f}", "phi": f"{phi:.1f}"}
        for isomer in ISOMERS:
            folder = runs_dir / f"{isomer}_th{int(theta):d}_ph{int(phi):d}"
            out_files = list(folder.glob(f"*_qed_ccsd_input.out"))
            if len(out_files) != 1:
                raise FileNotFoundError(f"{folder}: expected exactly one .out file, found {len(out_files)}")
            text = out_files[0].read_text()
            scf = SCF_RE.findall(text)
            ccsd = CCSD_RE.findall(text)
            if len(scf) != 1 or len(ccsd) != 1:
                raise ValueError(
                    f"{out_files[0]}: expected 1 SCF + 1 CCSD energy, found "
                    f"{len(scf)} SCF + {len(ccsd)} CCSD"
                )
            row[f"E_SCF_{isomer}_int"] = scf[0]
            row[f"E_CCSD_{isomer}_int"] = ccsd[0]
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build relaxed_qed_ccsd_intermediate_scans.csv from run outputs."
    )
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Folder containing the ortho/meta/para th* ph* run subfolders.",
    )
    parser.add_argument(
        "--grid",
        type=Path,
        default=Path(__file__).resolve().parent.parent.parent / "unrelaxed_qed_ccsd_intermediate_scans.csv",
        help="Reference CSV supplying the theta/phi grid (and output location).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV path (defaults to <grid parent>/relaxed_qed_ccsd_intermediate_scans.csv).",
    )
    args = parser.parse_args()

    grid_path = args.grid.resolve()
    if args.output is None:
        output_path = grid_path.parent / "relaxed_qed_ccsd_intermediate_scans.csv"
    else:
        output_path = args.output.resolve()

    rows = parse_runs_dir(args.runs_dir.resolve(), grid_path)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELD_NAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {output_path}")


if __name__ == "__main__":
    main()