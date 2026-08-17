#!/usr/bin/env python
"""
trajectory_campaign_no_freq.py

Optimization-only companion to trajectory_campaign.py for the direction-sampling
sweep. It reads the same trajectory_sampling_direction_*.dat files, validates the
same fixed physics parameters against cavity_common.py, and uses the same
campaign.optimize_cell restart loop -- but stops after the geometry relaxation.

This is meant for a faster, looser "relax the energy" pass where no Hessian, ZPE,
imaginary-mode check, or stalled-geometry promotion to a frequency calculation is
desired.

Usage:
    python trajectory_campaign_no_freq.py --isomer ortho
    python trajectory_campaign_no_freq.py --isomer ortho --dry-run
    python trajectory_campaign_no_freq.py --isomer ortho --list-ids
    python trajectory_campaign_no_freq.py --isomer ortho --summarize
    python trajectory_campaign_no_freq.py --isomer meta \
        --directions trajectory_sampling_direction_B.dat

Embarrassingly-parallel run:
    export PSI4_THREADS=4
    python trajectory_campaign_no_freq.py --isomer ortho --list-ids | \
        parallel -j 4 python trajectory_campaign_no_freq.py --isomer ortho --only {} --no-summarize
    python trajectory_campaign_no_freq.py --isomer ortho --summarize
"""

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

from campaign import (
    CONV_THRESHOLD,
    ISOMER_XYZ,
    optimize_cell,
    read_last_xyz_frame,
    save_json,
    _maybe_json,
)
from trajectory_campaign import (
    DEFAULT_DIRECTIONS_FILE,
    parse_directions_file,
    validate_fixed_params,
    direction_set_tag,
    build_manifest,
)


BASE_DIR = Path(__file__).resolve().parent
NO_FREQ_RUNS_DIR = BASE_DIR / "runs_trajectory_no_freq"


def run_cell(cell, force=False, reoptimize=False, continue_stalled=False):
    """Run one direction cell through optimization only."""
    import cavity_common as cc

    cell_dir = NO_FREQ_RUNS_DIR / cell["id"]
    cell_dir.mkdir(parents=True, exist_ok=True)

    done_marker = cell_dir / "DONE"
    if done_marker.exists() and not force:
        print(f"[skip] {cell['id']} (already DONE)")
        return "done"

    lambda_vector = cc.lambda_vector_for(cell["theta"], cell["phi"], cell["magnitude"])
    save_json(cell_dir / "cell.json", {**cell, "lambda_vector": lambda_vector})

    opt_status_path = cell_dir / "opt_status.json"
    opt_status = None
    if opt_status_path.exists() and not force:
        opt_status = json.loads(opt_status_path.read_text())

    already_converged = bool(opt_status and opt_status.get("converged"))
    continuable = bool(
        continue_stalled and opt_status and not already_converged and not reoptimize
        and (cell_dir / "opt_traj.xyz").exists()
    )

    if already_converged and not reoptimize:
        pass
    elif continuable:
        print(f"[continue] {cell['id']}: resuming from last trajectory geometry")
        symbols, coords_ang = read_last_xyz_frame(cell_dir / "opt_traj.xyz")
        opt_status = optimize_cell(
            cell, cell_dir, lambda_vector, resume_from=(symbols, coords_ang)
        )
    else:
        print(f"[opt ] {cell['id']} ...")
        opt_status = optimize_cell(cell, cell_dir, lambda_vector)

    if not opt_status.get("converged"):
        print(
            f"[STALL] {cell['id']}: not converged after "
            f"{opt_status.get('attempts')} attempts "
            f"(|g|={opt_status.get('final_gnorm')})."
        )
        return "stalled"

    done_marker.write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    print(
        f"[done] {cell['id']}: E = {opt_status.get('final_energy_hartree'):.10f} Ha, "
        f"|g| = {opt_status.get('final_gnorm'):.3e}"
    )
    return "done"


def summarize(cells, tag, isomer):
    """Scan runs_trajectory_no_freq/ and build an optimization-only rollup."""
    rows = []
    for cell in cells:
        cell_dir = NO_FREQ_RUNS_DIR / cell["id"]
        opt = _maybe_json(cell_dir / "opt_status.json")
        rows.append({
            "cell_id": cell["id"],
            "isomer": cell["isomer"],
            "direction_set": cell["direction_set"],
            "direction_row": cell["direction_row"],
            "theta": cell["theta"],
            "phi": cell["phi"],
            "magnitude": cell["magnitude"],
            "opt_converged": opt.get("converged") if opt else None,
            "attempts": opt.get("attempts") if opt else None,
            "final_gnorm": opt.get("final_gnorm") if opt else None,
            "E_opt_hartree": opt.get("final_energy_hartree") if opt else None,
            "wall_seconds": opt.get("wall_seconds") if opt else None,
        })

    NO_FREQ_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    base = f"trajectory_campaign_no_freq_{isomer}_dir{tag}_status"
    save_json(NO_FREQ_RUNS_DIR / f"{base}.json", rows)

    cols = [
        "cell_id", "isomer", "direction_set", "direction_row", "theta", "phi",
        "magnitude", "opt_converged", "attempts", "final_gnorm",
        "E_opt_hartree", "wall_seconds",
    ]
    csv_path = NO_FREQ_RUNS_DIR / f"{base}.csv"
    with open(csv_path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {NO_FREQ_RUNS_DIR / (base + '.json')}")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Direction-sweep optimization-only campaign driver."
    )
    parser.add_argument("--isomer", required=True, choices=sorted(ISOMER_XYZ),
                        help="isomer to run for every direction in the file")
    parser.add_argument("--directions", metavar="FILE",
                        default=DEFAULT_DIRECTIONS_FILE,
                        help=f"directions file (default: {DEFAULT_DIRECTIONS_FILE})")
    parser.add_argument("--dry-run", action="store_true", help="print the manifest and exit")
    parser.add_argument("--only", metavar="SUBSTR", help="only run cells whose id contains SUBSTR")
    parser.add_argument("--list-ids", action="store_true",
                        help="print one cell id per line (honors --only) and exit")
    parser.add_argument("--summarize", action="store_true",
                        help="rebuild the optimization rollup CSV/JSON and exit")
    parser.add_argument("--no-summarize", action="store_true",
                        help="skip the end-of-run rollup")
    parser.add_argument("--force", action="store_true",
                        help="redo cells even if already DONE")
    parser.add_argument("--reoptimize", action="store_true",
                        help="re-run the optimization from scratch")
    parser.add_argument("--continue-stalled", action="store_true",
                        help="resume a stalled optimization from the last opt_traj.xyz geometry")
    args = parser.parse_args()

    directions_path = Path(args.directions)
    if not directions_path.is_absolute():
        directions_path = BASE_DIR / directions_path

    rows, fixed = parse_directions_file(directions_path)
    tag = direction_set_tag(directions_path)
    all_cells = build_manifest(args.isomer, rows, directions_path)

    cells = all_cells
    if args.only:
        cells = [c for c in cells if args.only in c["id"]]

    if args.list_ids:
        for c in cells:
            print(c["id"])
        return

    if args.dry_run:
        print(f"Directions file: {directions_path.name} ({len(rows)} directions, set '{tag}')")
        print(
            f"Fixed: |lambda|={fixed['lambda_mag']} omega={fixed['omega']} "
            f"charge={fixed['charge']} mult={fixed['multiplicity']} "
            f"{fixed['functional']}/{fixed['basis']}"
        )
        print(f"{len(cells)} optimization-only cells:")
        for c in cells:
            print(
                f"  {c['id']:40s} theta={c['theta']:>6g} phi={c['phi']:>6g} "
                f"|lambda|={c['magnitude']:.2f}"
            )
        return

    validate_fixed_params(fixed, directions_path.name)

    if args.summarize:
        summarize(all_cells, tag, args.isomer)
        return

    NO_FREQ_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    print(
        f"Running {len(cells)} optimization-only cells "
        f"(isomer={args.isomer}, direction set '{tag}', |lambda|={fixed['lambda_mag']}). "
        f"Runs dir: {NO_FREQ_RUNS_DIR}"
    )

    counts = {"done": 0, "stalled": 0, "error": 0}
    for c in cells:
        try:
            outcome = run_cell(
                c, force=args.force, reoptimize=args.reoptimize,
                continue_stalled=args.continue_stalled,
            )
            counts[outcome] = counts.get(outcome, 0) + 1
        except Exception:
            counts["error"] += 1
            print(f"[ERROR] {c['id']}:\n{traceback.format_exc()}", file=sys.stderr)

    if not args.no_summarize and len(cells) > 1:
        summarize(all_cells, tag, args.isomer)
    print(f"\nOptimization-only campaign finished: {counts}")


if __name__ == "__main__":
    main()
