#!/usr/bin/env python
"""
trajectory_campaign.py

Companion driver to campaign.py for the DIRECTION-sampling sweep: the lambda
magnitude and all other physics parameters stay fixed, while the cavity
coupling direction (theta, phi) varies, read from a tab-separated file such as
trajectory_sampling_direction_A.dat:

    #Theta  Phi  lambda_mag  omega  charge  multiplicity  functional  basis
    101     31   0.1         0.06615  1     1             wb97x       6-311G*
    ...

Design decisions:
  * One isomer per invocation (--isomer ortho|meta|para); the manifest is
    (direction rows) x (that isomer).
  * Each cell runs the SAME optimize -> gate -> Hessian -> ZPE pipeline as
    campaign.py. The per-cell steps (optimize_cell, frequency_cell,
    promote_stalled) are imported from campaign.py, not re-implemented, so the
    convergence gate, restart loop, promotion logic, and provenance files stay
    identical and there is exactly one copy of that consequential code.
  * cavity_common.py remains the single source of truth for the physics
    parameters. The columns lambda_mag/omega/charge/multiplicity/functional/
    basis in the directions file are treated as a cross-check: every row must
    be identical, and the values must match cavity_common (and the fixed
    LAMBDA_MAG below), or the campaign refuses to run. The lambda vector itself
    is still built solely by cavity_common.lambda_vector_for.
  * Runs are archived under runs_trajectory/<cell_id>/ (separate from
    campaign.py's runs/) with the same DONE-marker resume semantics.

Usage:
    python trajectory_campaign.py --isomer ortho                # run/resume
    python trajectory_campaign.py --isomer ortho --dry-run      # print manifest
    python trajectory_campaign.py --isomer ortho --list-ids     # ids for a launcher
    python trajectory_campaign.py --isomer ortho --summarize    # rebuild rollup
    python trajectory_campaign.py --isomer meta \
        --directions trajectory_sampling_direction_B.dat        # another direction set
   Add --create-input-scripts flag to explicitly write input scripts for checking
Embarrassingly-parallel run (cap concurrency so cells don't oversubscribe):
    export PSI4_THREADS=4
    python trajectory_campaign.py --isomer ortho --list-ids | \
        parallel -j 4 python trajectory_campaign.py --isomer ortho --only {} --no-summarize
    python trajectory_campaign.py --isomer ortho --summarize
"""

import argparse
import json
import sys
import traceback
from pathlib import Path

# Reuse the consequential per-cell machinery from campaign.py verbatim.
from campaign import (
    ISOMER_XYZ,
    CONV_THRESHOLD,
    optimize_cell,
    frequency_cell,
    promote_stalled,
    save_json,
    _maybe_json,
    read_xyz,
    read_last_xyz_frame,
)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
TRAJ_RUNS_DIR = BASE_DIR / "runs_trajectory"

DEFAULT_DIRECTIONS_FILE = "trajectory_sampling_direction_A.dat"

# The fixed coupling magnitude for the whole direction sweep. Cross-checked
# against the lambda_mag column of the directions file at parse time.
LAMBDA_MAG = 0.10

# Expected column order in the directions file (after the '#' header line).
_COLUMNS = ("theta", "phi", "lambda_mag", "omega", "charge",
            "multiplicity", "functional", "basis")


# -----------------------------------------------------------------------------
# Directions file: parse + validate
# -----------------------------------------------------------------------------

def parse_directions_file(path):
    """
    Read a trajectory_sampling_direction_*.dat file.

    Returns (rows, fixed) where rows is a list of dicts with keys _COLUMNS plus
    'row' (1-based data-row index), and fixed is the single dict of fixed
    parameters shared by every row. Raises ValueError on any malformed line or
    if the fixed parameters are not identical across rows.
    """
    path = Path(path)
    rows = []
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) != len(_COLUMNS):
                raise ValueError(
                    f"{path.name}:{lineno}: expected {len(_COLUMNS)} columns "
                    f"({', '.join(_COLUMNS)}), got {len(parts)}: {stripped!r}"
                )
            try:
                row = {
                    "theta": float(parts[0]),
                    "phi": float(parts[1]),
                    "lambda_mag": float(parts[2]),
                    "omega": float(parts[3]),
                    "charge": int(parts[4]),
                    "multiplicity": int(parts[5]),
                    "functional": parts[6],
                    "basis": parts[7],
                }
            except ValueError as exc:
                raise ValueError(f"{path.name}:{lineno}: bad value ({exc})") from exc
            row["row"] = len(rows) + 1
            rows.append(row)

    if not rows:
        raise ValueError(f"{path.name}: no data rows found")

    # Fixed parameters must be identical on every row.
    fixed_keys = ("lambda_mag", "omega", "charge", "multiplicity",
                  "functional", "basis")
    fixed = {k: rows[0][k] for k in fixed_keys}
    for r in rows:
        for k in fixed_keys:
            if r[k] != fixed[k]:
                raise ValueError(
                    f"{path.name}: fixed parameter '{k}' varies across rows "
                    f"({fixed[k]!r} vs {r[k]!r} on data row {r['row']}). "
                    f"This campaign holds everything but (theta, phi) fixed."
                )
    return rows, fixed


def validate_fixed_params(fixed, directions_name):
    """
    Refuse to run if the directions file disagrees with cavity_common (the
    single source of truth) or with the campaign's fixed LAMBDA_MAG.
    """
    import cavity_common as cc

    problems = []

    def check(name, file_val, expected, eq):
        if not eq(file_val, expected):
            problems.append(f"  {name}: file has {file_val!r}, expected {expected!r}")

    check("lambda_mag", fixed["lambda_mag"], LAMBDA_MAG,
          lambda a, b: abs(a - b) < 1e-12)
    check("omega", fixed["omega"], cc.OMEGA, lambda a, b: abs(a - b) < 1e-12)
    check("charge", fixed["charge"], cc.CHARGE, lambda a, b: a == b)
    check("multiplicity", fixed["multiplicity"], cc.MULTIPLICITY,
          lambda a, b: a == b)
    check("functional", fixed["functional"], cc.FUNCTIONAL,
          lambda a, b: a.lower() == b.lower())
    check("basis", fixed["basis"], cc.PSI4_OPTIONS["basis"],
          lambda a, b: a.lower() == b.lower())

    if problems:
        raise SystemExit(
            f"ABORT: {directions_name} disagrees with cavity_common.py / "
            f"LAMBDA_MAG:\n" + "\n".join(problems) +
            "\nEither fix the file or update the single source of truth "
            "deliberately -- this campaign will not silently override it."
        )


# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------

def direction_set_tag(directions_path):
    """
    Short tag identifying the direction set, used in cell ids and rollup names.
    'trajectory_sampling_direction_A.dat' -> 'A'; otherwise the file stem.
    """
    stem = Path(directions_path).stem
    prefix = "trajectory_sampling_direction_"
    return stem[len(prefix):] if stem.startswith(prefix) and len(stem) > len(prefix) else stem


def _fmt_angle(x):
    """Format an angle for use in a cell id: 101.0 -> '101', 113.5 -> '113.5'."""
    return f"{x:g}"


def build_manifest(isomer, rows, directions_path):
    """Expand the direction rows into explicit cells for one isomer."""
    tag = direction_set_tag(directions_path)
    cells = []
    for r in rows:
        cells.append({
            # Row index in the id keeps ids unique even if a (theta, phi)
            # pair repeats in the file, and preserves file order under sort.
            "id": (f"{isomer}_dir{tag}_r{r['row']:03d}"
                   f"_th{_fmt_angle(r['theta'])}_ph{_fmt_angle(r['phi'])}"),
            "isomer": isomer,
            "xyz": ISOMER_XYZ[isomer],
            "theta": r["theta"],
            "phi": r["phi"],
            "magnitude": r["lambda_mag"],
            "direction_set": tag,
            "direction_row": r["row"],
        })
    return cells


# -----------------------------------------------------------------------------
# Stand-alone input-script generation (--create-input-scripts)
# -----------------------------------------------------------------------------

INPUT_SCRIPTS_DIR = BASE_DIR / "input_scripts"

# Stand-alone script template. Deliberately plain string substitution (not an
# f-string / .format()) so the many literal braces in the generated script's
# own dicts and f-strings don't need escaping. Every dynamic value is a
# unique __TOKEN__ replaced in build_input_script().
_INPUT_SCRIPT_TEMPLATE = '''import numpy as np
import psi4
from cqed_scf import CQEDCalculator
from cqed_scf.drivers import bfgs_optimize
from cqed_scf.utils import write_xyz, ANGSTROM_TO_BOHR, generate_field_vector_from_theta_and_phi
# =========================
# Psi4 geometry (angstrom)
# =========================
geometry_string = """
__GEOMETRY__"""
# =========================
# Psi4 options
# =========================
psi4_options = {
    "basis": "__BASIS__",
    "scf_type": "df",
    "e_convergence": 1e-12,
    "d_convergence": 1e-12,
    "dft_radial_points": 99,
    "dft_spherical_points": 590,
    "dft_pruning_scheme": "none"
}
psi4.set_options(psi4_options)
# =========================
# Cavity parameters
# =========================
theta = __THETA__  # degrees from z-axis
phi = __PHI__  # degrees from x-axis in xy-plane
field_vector_center = generate_field_vector_from_theta_and_phi(theta, phi)
print(F"Field Magnitude before scaling is {np.linalg.norm(field_vector_center)}")
print(F"Field Vector is")
print(field_vector_center)
lambda_direction = np.asarray(field_vector_center, dtype=float)
lambda_direction /= np.linalg.norm(lambda_direction)
lam_mag = __LAM_MAG__
omega = __OMEGA__
mol = psi4.geometry(geometry_string)
symbols = [mol.symbol(i) for i in range(mol.natom())]
lambda_vector = (lam_mag * lambda_direction).tolist()
xyz_file = "__XYZ_FILE__"
print(F"Mag after scaling {np.linalg.norm(lambda_vector)}")
open(xyz_file, "w").close()
print(f"Running optimization for |lambda| = {lam_mag:.2f} with vector {lambda_vector}")
calc = CQEDCalculator(
    lambda_vector=lambda_vector,
    psi4_options=psi4_options,
    omega=omega,
    density_fitting=True,
    charge=__CHARGE__,
    multiplicity=__MULTIPLICITY__,
    functional="__FUNCTIONAL__",
)
opt_result, _ = bfgs_optimize(
    calculator=calc,
    geometry=geometry_string,
    canonical="psi4",
    gtol=1e-6,
    maxiter=50,
    debug=True,
    project_tr_rot=True,
    projection_debug=True,
)
coords_opt_angstrom = opt_result.x.reshape(-1, 3) / ANGSTROM_TO_BOHR
write_xyz(
    xyz_file,
    symbols,
    coords_opt_angstrom,
    comment=f"FINAL OPTIMIZED | lambda = {lam_mag:.2f} | E = {opt_result.fun:.10f} Ha",
    mode="a",
)
# =========================
# Summary
# =========================
print("\\nOptimization finished.")
print(f"Converged: {opt_result.success}")
print(f"Final energy (Ha): {opt_result.fun:.10f}")
print(f"XYZ trajectory written to: {xyz_file}")
'''


def build_input_script(cell, fixed):
    """
    Render the stand-alone input script for one cell as plain text.

    Pure string formatting -- no psi4/cqed_scf computation is performed here.
    The only repo import used is cavity_common.psi4_geometry_string (via
    campaign.read_xyz) to embed the geometry literally, per the single-source-
    of-truth policy for that formatting.
    """
    import cavity_common as cc

    symbols, coords0 = read_xyz(BASE_DIR / cell["xyz"])
    geometry_string = cc.psi4_geometry_string(symbols, coords0)
    xyz_file = f"opt_{cell['id']}.xyz"

    return (_INPUT_SCRIPT_TEMPLATE
            .replace("__GEOMETRY__", geometry_string)
            .replace("__BASIS__", fixed["basis"])
            .replace("__THETA__", _fmt_angle(cell["theta"]))
            .replace("__PHI__", _fmt_angle(cell["phi"]))
            .replace("__LAM_MAG__", repr(fixed["lambda_mag"]))
            .replace("__OMEGA__", repr(fixed["omega"]))
            .replace("__XYZ_FILE__", xyz_file)
            .replace("__CHARGE__", str(fixed["charge"]))
            .replace("__MULTIPLICITY__", str(fixed["multiplicity"]))
            .replace("__FUNCTIONAL__", fixed["functional"]))


def create_input_scripts(cells, fixed):
    """Write one stand-alone input script per cell into input_scripts/."""
    INPUT_SCRIPTS_DIR.mkdir(parents=True, exist_ok=True)
    for cell in cells:
        script = build_input_script(cell, fixed)
        out_path = INPUT_SCRIPTS_DIR / f"input_{cell['id']}.py"
        out_path.write_text(script)
        print(f"Wrote {out_path}")
    print(f"\n{len(cells)} input script(s) written to {INPUT_SCRIPTS_DIR}")


# -----------------------------------------------------------------------------
# Per-cell driver (mirrors campaign.run_cell, but archives under runs_trajectory/)
# -----------------------------------------------------------------------------

def run_cell(cell, force=False, reoptimize=False, continue_stalled=False):
    """Run one cell end-to-end with the same resume/checkpoint logic as campaign.py."""
    import cavity_common as cc

    cell_dir = TRAJ_RUNS_DIR / cell["id"]
    cell_dir.mkdir(parents=True, exist_ok=True)

    done_marker = cell_dir / "DONE"
    if done_marker.exists() and not force:
        print(f"[skip] {cell['id']} (already DONE)")
        return "done"

    # The ONLY place lambda is constructed: cavity_common, from (theta, phi, |lambda|).
    lambda_vector = cc.lambda_vector_for(cell["theta"], cell["phi"], cell["magnitude"])
    save_json(cell_dir / "cell.json", {**cell, "lambda_vector": lambda_vector})

    # ---- Optimization (identical gate/promotion logic to campaign.py) ----
    opt_status_path = cell_dir / "opt_status.json"
    opt_status = None
    if opt_status_path.exists() and not force:
        opt_status = json.loads(opt_status_path.read_text())

    already_converged = bool(opt_status and opt_status.get("converged"))
    promotable = bool(
        opt_status and not already_converged and not reoptimize
        and opt_status.get("final_gnorm") is not None
        and opt_status["final_gnorm"] < CONV_THRESHOLD
        and (cell_dir / "opt_traj.xyz").exists()
    )
    # A genuine stall: exhausted MAX_RESTARTS in campaign.optimize_cell without
    # ever getting under CONV_THRESHOLD (so it's not "promotable"). With
    # --continue-stalled, resume from the last trajectory geometry -- appending
    # fresh restart attempts -- instead of starting over from the original
    # isomer .xyz.
    continuable = bool(
        continue_stalled and opt_status and not already_converged and not promotable
        and (cell_dir / "opt_traj.xyz").exists()
    )

    if already_converged:
        pass
    elif promotable:
        print(f"[promote] {cell['id']}: |g|={opt_status['final_gnorm']:.2e} "
              f"< {CONV_THRESHOLD:.0e}; accepting stalled geometry, no re-opt")
        opt_status = promote_stalled(cell, cell_dir, opt_status)
    elif continuable:
        print(f"[continue] {cell['id']}: stalled at |g_proj|={opt_status['final_gnorm']:.2e} "
              f">= {CONV_THRESHOLD:.0e} after {opt_status.get('attempts')} attempts; "
              f"resuming from last trajectory geometry instead of the original isomer geometry")
        symbols, coords_ang = read_last_xyz_frame(cell_dir / "opt_traj.xyz")
        restart_geometry_string = cc.psi4_geometry_string(symbols, coords_ang)
        print(
            f"restarting stalled geometry optimization with direction "
            f"{cell['theta']},{cell['phi']} which had Energy E = "
            f"{opt_status.get('final_energy_hartree')} and |g_proj| = "
            f"{opt_status.get('final_gnorm')} with the following geometry:\n"
            f"{restart_geometry_string}"
        )
        opt_status = optimize_cell(cell, cell_dir, lambda_vector,
                                    resume_from=(symbols, coords_ang))
    else:
        print(f"[opt ] {cell['id']} ...")
        opt_status = optimize_cell(cell, cell_dir, lambda_vector)

    if not opt_status.get("converged"):
        print(f"[STALL] {cell['id']}: not converged after "
              f"{opt_status.get('attempts')} attempts "
              f"(|g|={opt_status.get('final_gnorm')}). Skipping frequencies.")
        return "stalled"

    # ---- Frequencies ----
    freq_path = cell_dir / "frequencies.json"
    if freq_path.exists() and not force:
        freq = json.loads(freq_path.read_text())
    else:
        print(f"[freq] {cell['id']} ...")
        freq = frequency_cell(cell, cell_dir, lambda_vector)

    import time
    done_marker.write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))

    flag = " <-- IMAGINARY MODE(S)" if freq.get("n_imaginary", 0) > 0 else ""
    print(f"[done] {cell['id']}: ZPE = {freq.get('zpe_hartree'):.8f} Ha, "
          f"n_imag = {freq.get('n_imaginary')}{flag}")
    return "done"


# -----------------------------------------------------------------------------
# Rollup / summary
# -----------------------------------------------------------------------------

def summarize(cells, tag, isomer):
    """Scan runs_trajectory/ and build a CSV + JSON rollup for this manifest."""
    rows = []
    for cell in cells:
        cell_dir = TRAJ_RUNS_DIR / cell["id"]
        opt = _maybe_json(cell_dir / "opt_status.json")
        freq = _maybe_json(cell_dir / "frequencies.json")
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
            "zpe_hartree": freq.get("zpe_hartree") if freq else None,
            "n_imaginary": freq.get("n_imaginary") if freq else None,
            "lowest_freq_cm": freq.get("lowest_freq_cm") if freq else None,
        })

    TRAJ_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    base = f"trajectory_campaign_{isomer}_dir{tag}_status"
    save_json(TRAJ_RUNS_DIR / f"{base}.json", rows)

    cols = ["cell_id", "isomer", "direction_set", "direction_row", "theta",
            "phi", "magnitude", "opt_converged", "attempts", "final_gnorm",
            "E_opt_hartree", "zpe_hartree", "n_imaginary", "lowest_freq_cm"]
    csv_path = TRAJ_RUNS_DIR / f"{base}.csv"
    with open(csv_path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {TRAJ_RUNS_DIR / (base + '.json')}")
    return rows


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Direction-sweep (fixed |lambda|) opt + frequency campaign driver.")
    parser.add_argument("--isomer", required=True, choices=sorted(ISOMER_XYZ),
                        help="isomer to run for every direction in the file")
    parser.add_argument("--directions", metavar="FILE",
                        default=DEFAULT_DIRECTIONS_FILE,
                        help=f"directions file (default: {DEFAULT_DIRECTIONS_FILE})")
    parser.add_argument("--dry-run", action="store_true", help="print the manifest and exit")
    parser.add_argument("--only", metavar="SUBSTR", help="only run cells whose id contains SUBSTR")
    parser.add_argument("--list-ids", action="store_true",
                        help="print one cell id per line (honors --only) and exit")
    parser.add_argument("--create-input-scripts", action="store_true",
                        help="write one stand-alone Python input script per direction "
                             "row into input_scripts/ (honors --only) and exit")
    parser.add_argument("--summarize", action="store_true", help="rebuild the rollup CSV/JSON and exit")
    parser.add_argument("--no-summarize", action="store_true",
                        help="skip the end-of-run rollup (for parallel single-cell shards)")
    parser.add_argument("--force", action="store_true", help="redo cells even if already DONE")
    parser.add_argument("--reoptimize", action="store_true",
                        help="re-run the optimization instead of promoting a stalled "
                             "geometry that already clears CONV_THRESHOLD")
    parser.add_argument("--continue-stalled", action="store_true",
                        help="for a cell that exhausted MAX_RESTARTS without converging, "
                             "resume optimization from the last opt_traj.xyz geometry "
                             "(appending further restart attempts) instead of starting "
                             "over from the original isomer geometry")
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
        print(f"Directions file: {directions_path.name} ({len(rows)} directions, "
              f"set '{tag}')")
        print(f"Fixed: |lambda|={fixed['lambda_mag']} omega={fixed['omega']} "
              f"charge={fixed['charge']} mult={fixed['multiplicity']} "
              f"{fixed['functional']}/{fixed['basis']}")
        print(f"{len(cells)} cells:")
        for c in cells:
            print(f"  {c['id']:40s} theta={c['theta']:>6g} phi={c['phi']:>6g} "
                  f"|lambda|={c['magnitude']:.2f}")
        return

    # Physics gate: file must agree with cavity_common before any compute.
    validate_fixed_params(fixed, directions_path.name)

    if args.create_input_scripts:
        create_input_scripts(cells, fixed)
        return

    if args.summarize:
        summarize(all_cells, tag, args.isomer)
        return

    TRAJ_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Running {len(cells)} cells (isomer={args.isomer}, direction set "
          f"'{tag}', |lambda|={fixed['lambda_mag']}). Runs dir: {TRAJ_RUNS_DIR}")

    counts = {"done": 0, "stalled": 0, "error": 0}
    for c in cells:
        try:
            outcome = run_cell(c, force=args.force, reoptimize=args.reoptimize,
                               continue_stalled=args.continue_stalled)
            counts[outcome] = counts.get(outcome, 0) + 1
        except Exception:
            counts["error"] += 1
            print(f"[ERROR] {c['id']}:\n{traceback.format_exc()}", file=sys.stderr)
            # Keep going -- one bad cell shouldn't sink a multi-day campaign.

    if not args.no_summarize and len(cells) > 1:
        summarize(all_cells, tag, args.isomer)
    print(f"\nCampaign finished: {counts}")


if __name__ == "__main__":
    main()

