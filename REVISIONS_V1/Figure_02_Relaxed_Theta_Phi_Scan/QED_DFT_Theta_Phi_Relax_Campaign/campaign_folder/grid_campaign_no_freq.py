#!/usr/bin/env python
"""
grid_campaign_no_freq.py

Constrained QED-DFT geometry-relaxation campaign over the (theta, phi) grid
in intermediate_scans.csv (the same grid used for the QED-CCSD/EOM-CCSD
theta/phi energy scan), for one Wheland-intermediate isomer at a time.

Unlike trajectory_campaign.py / trajectory_campaign_no_freq.py (which treat
each (theta, phi) direction as an independent cell, safe to run in any order
or in parallel), this campaign chains a WARM START across phi within each
theta ring: for a fixed theta, phi=0 is optimized starting from the isomer's
pristine geometry; phi=<next> is then optimized starting from phi=0's
converged geometry; and so on through every phi listed for that theta in
intermediate_scans.csv. Each theta ring is an independent chain -- it always
starts fresh from the pristine isomer.xyz at phi=0, never from a neighboring
theta ring -- and different theta rings (and different isomers) are fully
independent of each other and safe to run concurrently.

Differences from trajectory_campaign_no_freq.py:
  * Directions come from intermediate_scans.csv's (theta, phi) grid, not a
    trajectory_sampling_direction_*.dat file. Only the (theta, phi) columns
    are used -- the energy columns in that CSV are a different level of
    theory's precomputed scan and are not consulted here.
  * The convergence gate is tighter: |g_projected| < GRID_CONV_THRESHOLD
    (2e-3), passed explicitly to campaign.optimize_cell's new conv_threshold
    parameter -- this campaign does not change the default used by
    campaign.py / trajectory_campaign.py.
  * No frequency/Hessian step at all (same as trajectory_campaign_no_freq.py).
  * On every converged cell, also writes a
    <cell_id>_qed_ccsd_input.json (e.g. ortho_th63_ph126_qed_ccsd_input.json;
    via ccsd_input.py) so the relaxed geometry can be handed straight to a
    downstream QED-CCSD calculation at the same (theta, phi) direction and a
    fixed lambda magnitude of 0.1.
  * If a phi step in a ring fails to converge within MAX_RESTARTS, that ring
    HALTS at that point -- later phi values in the same ring are left
    pending rather than warm-started from a geometry that isn't actually
    converged. Re-running later (after investigating/fixing the stall)
    resumes the ring from wherever it left off.

The unit of parallelism is the (isomer, theta) RING, not the individual
cell -- phi steps within one ring are strictly sequential.

Usage:
    python grid_campaign_no_freq.py --isomer ortho                # run/resume all rings
    python grid_campaign_no_freq.py --isomer ortho --dry-run       # print the manifest
    python grid_campaign_no_freq.py --isomer ortho --list-groups   # ring ids, for a launcher
    python grid_campaign_no_freq.py --isomer ortho --only-group th63
    python grid_campaign_no_freq.py --isomer ortho --summarize

Embarrassingly-parallel run (cap concurrency so cells don't oversubscribe):
    export PSI4_THREADS=4
    python grid_campaign_no_freq.py --isomer ortho --list-groups | \\
        parallel -j 4 python grid_campaign_no_freq.py --isomer ortho --only-group {} --no-summarize
    python grid_campaign_no_freq.py --isomer ortho --summarize
"""

import argparse
import csv
import sys
import time
import traceback
from pathlib import Path

from campaign import (
    ISOMER_XYZ,
    optimize_cell,
    read_xyz,
    save_json,
    _maybe_json,
)

import ccsd_input

BASE_DIR = Path(__file__).resolve().parent
GRID_RUNS_DIR = BASE_DIR / "runs_grid_no_freq"

DEFAULT_SCAN_FILE = "intermediate_scans.csv"

# Tighter than campaign.py's CONV_THRESHOLD (2.5e-3) -- a deliberately loose
# "light relaxation" gate for this sweep, per Jay's spec.
GRID_CONV_THRESHOLD = 2e-3

# Fixed coupling magnitude for every cell in this campaign (matches the
# qed_lambdas value written into every generated QED-CCSD input file).
LAMBDA_MAG = 0.1


# -----------------------------------------------------------------------------
# Scan file: parse (theta, phi) grid
# -----------------------------------------------------------------------------

def parse_scan_file(path):
    """
    Read intermediate_scans.csv and return {theta: [phi, phi, ...]} with theta
    keys ascending and each theta's phi list ascending (starting at 0), which
    is the exact warm-start chain order. Only the theta,phi columns are used.
    """
    path = Path(path)
    by_theta = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if "theta" not in reader.fieldnames or "phi" not in reader.fieldnames:
            raise ValueError(
                f"{path.name}: expected 'theta' and 'phi' columns, got {reader.fieldnames}"
            )
        for row in reader:
            theta = float(row["theta"])
            phi = float(row["phi"])
            by_theta.setdefault(theta, []).append(phi)

    if not by_theta:
        raise ValueError(f"{path.name}: no data rows found")

    return {theta: sorted(set(phis)) for theta, phis in sorted(by_theta.items())}


def _fmt_angle(x):
    """Format an angle for use in an id: 63.0 -> '63', 31.5 -> '31.5'."""
    return f"{x:g}"


# -----------------------------------------------------------------------------
# Manifest: (isomer, theta) rings, each an ordered phi chain
# -----------------------------------------------------------------------------

def build_groups(isomer, by_theta):
    """
    Expand {theta: [phi...]} into an ordered list of ring dicts for one
    isomer. Each ring is {"id", "isomer", "theta", "cells": [cell, ...]},
    with cells in warm-start (phi-ascending) order.
    """
    groups = []
    for theta in sorted(by_theta):
        phis = by_theta[theta]
        group_id = f"{isomer}_th{_fmt_angle(theta)}"
        cells = []
        for i, phi in enumerate(phis):
            cells.append({
                "id": f"{isomer}_th{_fmt_angle(theta)}_ph{_fmt_angle(phi)}",
                "isomer": isomer,
                "xyz": ISOMER_XYZ[isomer],
                "theta": theta,
                "phi": phi,
                "magnitude": LAMBDA_MAG,
                "group_id": group_id,
                "chain_index": i,
                "is_first_in_ring": i == 0,
            })
        groups.append({"id": group_id, "isomer": isomer, "theta": theta, "cells": cells})
    return groups


# -----------------------------------------------------------------------------
# Per-cell / per-ring execution
# -----------------------------------------------------------------------------

def _cell_dir(cell):
    return GRID_RUNS_DIR / cell["id"]


def run_cell(cell, resume_from, force=False):
    """
    Run (or skip, if already DONE) one cell. Returns
    (outcome, symbols, coords_angstrom) where outcome is one of
    "done" (already done, loaded from disk), "converged" (just optimized and
    converged), "stalled" (just optimized, did not converge), or "error".
    symbols/coords_angstrom are the geometry to seed the NEXT phi in the ring
    (None if the cell did not converge).
    """
    import cavity_common as cc

    cell_dir = _cell_dir(cell)
    cell_dir.mkdir(parents=True, exist_ok=True)

    done_marker = cell_dir / "DONE"
    if done_marker.exists() and not force:
        symbols, coords_ang = read_xyz(cell_dir / "optimized.xyz")
        return "done", symbols, coords_ang

    lambda_vector = cc.lambda_vector_for(cell["theta"], cell["phi"], cell["magnitude"])
    save_json(cell_dir / "cell.json", {**cell, "lambda_vector": lambda_vector})

    print(f"[opt ] {cell['id']} ...")
    opt_status = optimize_cell(
        cell, cell_dir, lambda_vector,
        resume_from=resume_from,
        conv_threshold=GRID_CONV_THRESHOLD,
    )

    if not opt_status.get("converged"):
        print(f"[STALL] {cell['id']}: not converged after "
              f"{opt_status.get('attempts')} attempts "
              f"(|g|={opt_status.get('final_gnorm')}). Halting this ring.")
        return "stalled", None, None

    symbols, coords_ang = read_xyz(cell_dir / "optimized.xyz")

    # Named after the cell (<isomer>_th<theta>_ph<phi>), not just
    # "qed_ccsd_input.json", so a human moving these files around on a
    # different machine later (out of the runs_grid_no_freq/<cell_id>/
    # directory structure) can still tell which isomer/direction each one is.
    ccsd_path = cell_dir / f"{cell['id']}_qed_ccsd_input.json"
    ccsd_input.write_qed_ccsd_input(ccsd_path, cell["theta"], cell["phi"], symbols, coords_ang)

    done_marker.write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    print(f"[done] {cell['id']}: E = {opt_status.get('final_energy_hartree'):.10f} Ha, "
          f"|g| = {opt_status.get('final_gnorm'):.3e}, wrote {ccsd_path.name}")
    return "converged", symbols, coords_ang


def run_group(group, force=False, only=None):
    """
    Run one (isomer, theta) ring end-to-end: phi=0 seeded from the pristine
    isomer geometry, each subsequent phi seeded from the previous phi's
    converged geometry. Halts (leaves remaining cells pending) at the first
    stall. Returns a dict of per-outcome counts for this ring.
    """
    counts = {"done": 0, "converged": 0, "stalled": 0, "skipped": 0, "error": 0}
    prev_geom = None  # (symbols, coords_angstrom) of the last converged/loaded cell

    for cell in group["cells"]:
        if only and only not in cell["id"]:
            counts["skipped"] += 1
            # Still need this cell's geometry to seed the next one, if present.
            cell_dir = _cell_dir(cell)
            if (cell_dir / "DONE").exists():
                prev_geom = read_xyz(cell_dir / "optimized.xyz")
            continue

        resume_from = None if cell["is_first_in_ring"] else prev_geom
        if not cell["is_first_in_ring"] and prev_geom is None:
            print(f"[HALT ] {group['id']}: predecessor of {cell['id']} never converged; "
                  f"leaving remaining phi in this ring pending.")
            break

        try:
            outcome, symbols, coords_ang = run_cell(cell, resume_from, force=force)
        except Exception:
            counts["error"] += 1
            print(f"[ERROR] {cell['id']}:\n{traceback.format_exc()}", file=sys.stderr)
            break  # don't chain a warm start off of a cell that errored

        counts[outcome if outcome in counts else "error"] = counts.get(outcome, 0) + 1
        if outcome in ("done", "converged"):
            prev_geom = (symbols, coords_ang)
        else:
            # stalled: no good geometry to hand to the next phi -- halt the ring.
            break

    return counts


# -----------------------------------------------------------------------------
# Rollup / summary
# -----------------------------------------------------------------------------

def summarize(groups, isomer):
    """Scan runs_grid_no_freq/ and build a CSV + JSON rollup for this isomer."""
    rows = []
    for group in groups:
        for cell in group["cells"]:
            cell_dir = _cell_dir(cell)
            opt = _maybe_json(cell_dir / "opt_status.json")
            has_ccsd_input = (cell_dir / f"{cell['id']}_qed_ccsd_input.json").exists()
            rows.append({
                "cell_id": cell["id"],
                "group_id": cell["group_id"],
                "isomer": cell["isomer"],
                "theta": cell["theta"],
                "phi": cell["phi"],
                "chain_index": cell["chain_index"],
                "opt_converged": opt.get("converged") if opt else None,
                "attempts": opt.get("attempts") if opt else None,
                "final_gnorm": opt.get("final_gnorm") if opt else None,
                "E_opt_hartree": opt.get("final_energy_hartree") if opt else None,
                "conv_threshold": opt.get("conv_threshold") if opt else None,
                "qed_ccsd_input_written": has_ccsd_input,
            })

    GRID_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    base = f"grid_campaign_no_freq_{isomer}_status"
    save_json(GRID_RUNS_DIR / f"{base}.json", rows)

    cols = ["cell_id", "group_id", "isomer", "theta", "phi", "chain_index",
            "opt_converged", "attempts", "final_gnorm", "E_opt_hartree",
            "conv_threshold", "qed_ccsd_input_written"]
    csv_path = GRID_RUNS_DIR / f"{base}.csv"
    with open(csv_path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {GRID_RUNS_DIR / (base + '.json')}")
    return rows


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Theta/phi-grid constrained QED-DFT relaxation campaign "
                    "(no frequencies), warm-started across phi within each theta ring.")
    parser.add_argument("--isomer", required=True, choices=sorted(ISOMER_XYZ),
                        help="isomer to run for every (theta, phi) in the scan file")
    parser.add_argument("--scan-file", metavar="FILE", default=DEFAULT_SCAN_FILE,
                        help=f"CSV with theta,phi columns (default: {DEFAULT_SCAN_FILE})")
    parser.add_argument("--dry-run", action="store_true", help="print the manifest and exit")
    parser.add_argument("--list-groups", action="store_true",
                        help="print one ring (group) id per line (honors --only-group) and exit; "
                             "for piping to a job launcher -- each ring must run in ONE process, "
                             "since phi steps within it are sequential")
    parser.add_argument("--only-group", metavar="SUBSTR",
                        help="only run rings whose group id contains SUBSTR (e.g. 'th63')")
    parser.add_argument("--only", metavar="SUBSTR",
                        help="within the selected ring(s), only (re)run cells whose id contains "
                             "SUBSTR; other cells are used only as warm-start sources if already done")
    parser.add_argument("--force", action="store_true",
                        help="redo cells even if already DONE")
    parser.add_argument("--summarize", action="store_true", help="rebuild the rollup CSV/JSON and exit")
    parser.add_argument("--no-summarize", action="store_true",
                        help="skip the end-of-run rollup (use for parallel single-ring shards)")
    args = parser.parse_args()

    scan_path = Path(args.scan_file)
    if not scan_path.is_absolute():
        scan_path = BASE_DIR / scan_path

    by_theta = parse_scan_file(scan_path)
    all_groups = build_groups(args.isomer, by_theta)

    groups = all_groups
    if args.only_group:
        groups = [g for g in groups if args.only_group in g["id"]]

    if args.list_groups:
        for g in groups:
            print(g["id"])
        return

    if args.dry_run:
        n_cells = sum(len(g["cells"]) for g in all_groups)
        print(f"Scan file: {scan_path.name} ({len(all_groups)} theta rings, {n_cells} total cells)")
        print(f"Fixed: |lambda|={LAMBDA_MAG} conv_threshold={GRID_CONV_THRESHOLD:.1e} "
              f"isomer={args.isomer}")
        for g in groups:
            phis = ", ".join(_fmt_angle(c["phi"]) for c in g["cells"])
            print(f"  {g['id']:16s} theta={g['theta']:>5g}  {len(g['cells'])} phi steps: {phis}")
        return

    if args.summarize:
        summarize(all_groups, args.isomer)
        return

    GRID_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    n_cells = sum(len(g["cells"]) for g in groups)
    print(f"Running {len(groups)} theta ring(s), {n_cells} cells "
          f"(isomer={args.isomer}). Runs dir: {GRID_RUNS_DIR}")

    totals = {"done": 0, "converged": 0, "stalled": 0, "skipped": 0, "error": 0}
    for g in groups:
        counts = run_group(g, force=args.force, only=args.only)
        for k, v in counts.items():
            totals[k] = totals.get(k, 0) + v

    if not args.no_summarize and len(groups) > 1:
        summarize(all_groups, args.isomer)
    print(f"\nCampaign finished: {totals}")


if __name__ == "__main__":
    main()
