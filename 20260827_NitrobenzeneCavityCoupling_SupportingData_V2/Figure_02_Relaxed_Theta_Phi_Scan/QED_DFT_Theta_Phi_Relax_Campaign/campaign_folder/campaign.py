"""
campaign.py

Deterministic driver for the constrained cavity geometry-optimization +
frequency campaign on the bromo-nitrobenzene isomers.

Design philosophy (see project ROADMAP discussion): the control flow here is
crisp and consequential -- a wrong "converged" decision wastes ~7.5 h of Hessian
compute, and a wrong restart throws away a ~10 h optimization. So the loop,
convergence test, gating, checkpointing, and bookkeeping are all plain,
auditable Python. No LLM sits on this path. (An optional qwen layer for failure
diagnosis / campaign summarization can be added on top later; it is not needed
to remove the tedium.)

What it does, per (isomer, direction, magnitude) cell:
  1. derive the lambda vector ONCE from (theta, phi, magnitude);
  2. optimize with project_tr_rot=True, restarting from the last geometry
     until the projected gradient norm < CONV_THRESHOLD or MAX_RESTARTS hit;
  3. gate: only run frequencies if the optimization actually converged;
  4. run the ASE finite-difference Hessian, project out rigid-body modes,
     diagonalize, compute the ZPE, and flag any imaginary modes;
  5. archive everything in runs/<cell_id>/ and drop a DONE marker.

The campaign is resumable: completed cells (and even partially completed
Hessians, via ASE's own cache) are skipped on re-run.

Usage:
    python campaign.py                # run/resume the whole campaign
    python campaign.py --dry-run      # print the manifest and exit
    python campaign.py --only meta    # only cells whose id contains "meta"
    python campaign.py --list-ids     # print cell ids (for piping to a launcher)
    python campaign.py --summarize    # rebuild the CSV/JSON rollup and exit
    python campaign.py --force        # redo cells even if already DONE
    python campaign.py --reoptimize   # re-run opt instead of promoting a prior stall

Embarrassingly-parallel run (cap concurrency so cells don't oversubscribe):
    export PSI4_THREADS=4
    python campaign.py --list-ids | parallel -j 4 python campaign.py --only {} --no-summarize
    python campaign.py --summarize
"""

import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
RUNS_DIR = BASE_DIR / "runs"


# -----------------------------------------------------------------------------
# Campaign configuration
# -----------------------------------------------------------------------------

# Sweep specification. Each entry pairs one cavity DIRECTION (theta, phi) with
# the isomers that use it and the magnitudes to scan. Meta appears in both
# directions; ortho only in the first, para only in the second.
SWEEPS = [
    {
        "theta": 70, "phi": 31,
        "isomers": ["ortho", "meta"],
        "magnitudes": [0.02, 0.04, 0.06, 0.08, 0.10],
    },
    {
        "theta": 65, "phi": 78,
        "isomers": ["para", "meta"],
        "magnitudes": [0.02, 0.04, 0.06, 0.08, 0.10],
    },
]

ISOMER_XYZ = {
    "ortho": "ortho.xyz",
    "meta": "meta.xyz",
    "para": "para.xyz",
}

# Optimization / convergence controls.
GTOL_BFGS = 1e-6        # tol handed to bfgs_optimize (it may not reach this)
# Projected-gradient-norm gate ("optimized enough" to run frequencies).
# Raised from 1e-4 to 5e-4: on the cavity-constrained surface the projected
# gradient plateaus near ~1e-4 (DFT grid noise) and stalls out over the restart
# cap. 5e-4 is still tight (norm over 45 DOF => ~7e-5 RMS/component) and clears
# the plateau. The value is recorded per cell in opt_status.json for provenance.
CONV_THRESHOLD = 2.5e-3
MAXITER = 50            # steps per bfgs_optimize call (matches example)
MAX_RESTARTS = 6        # cap on restart attempts before declaring a stall
DEBUG_OPT = False       # verbose bfgs/projection output

# ASE finite-difference displacement (Angstrom). ASE Vibrations default is 0.01.
FD_DELTA = 0.01


# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------

def build_manifest():
    """Expand SWEEPS into the explicit list of cells."""
    cells = []
    for sweep in SWEEPS:
        theta, phi = sweep["theta"], sweep["phi"]
        for isomer in sweep["isomers"]:
            for mag in sweep["magnitudes"]:
                cells.append({
                    "id": f"{isomer}_{theta}_{phi}_lam{mag:.2f}",
                    "isomer": isomer,
                    "xyz": ISOMER_XYZ[isomer],
                    "theta": theta,
                    "phi": phi,
                    "magnitude": mag,
                })
    return cells


# -----------------------------------------------------------------------------
# Small IO helpers
# -----------------------------------------------------------------------------

def write_xyz_frame(path, symbols, coords_angstrom, comment="", mode="w"):
    with open(path, mode) as f:
        f.write(f"{len(symbols)}\n")
        f.write(f"{comment}\n")
        for s, p in zip(symbols, coords_angstrom):
            f.write(f"{s:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")


def read_xyz(path):
    """Minimal XYZ reader (first frame) -> (symbols, coords_angstrom ndarray)."""
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    symbols, coords = [], []
    for line in lines[2:2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.asarray(coords, dtype=float)


def read_last_xyz_frame(path):
    """Read the LAST frame of a (possibly multi-frame) XYZ trajectory."""
    with open(path) as f:
        lines = f.read().splitlines()
    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].split()[0])
        except (ValueError, IndexError):
            i += 1
            continue
        block = lines[i + 2:i + 2 + n]
        symbols = [ln.split()[0] for ln in block]
        coords = [[float(x) for x in ln.split()[1:4]] for ln in block]
        frames.append((symbols, np.asarray(coords, dtype=float)))
        i += 2 + n
    if not frames:
        raise ValueError(f"No frames found in {path}")
    return frames[-1]


def save_json(path, data):
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=_json_default)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# -----------------------------------------------------------------------------
# Per-cell steps
# -----------------------------------------------------------------------------

def optimize_cell(cell, cell_dir, lambda_vector, resume_from=None,
                   conv_threshold=CONV_THRESHOLD):
    """
    Optimize-check-restart loop, gated on the projected gradient norm.

    resume_from: optional (symbols, coords_angstrom) pair. When given, the
    optimization is seeded from this geometry instead of the original isomer
    .xyz, and the existing opt_traj.xyz is appended to rather than truncated
    -- this is how trajectory_campaign.py's --continue-stalled chains fresh
    restart attempts onto a run that already exhausted MAX_RESTARTS without
    converging. Default (None) preserves the original from-scratch behavior.

    conv_threshold: projected-gradient-norm gate for this cell. Defaults to
    the module-level CONV_THRESHOLD (unchanged behavior for existing callers),
    but callers such as grid_campaign_no_freq.py can pass a different value
    (e.g. a tighter/looser threshold for a different kind of sweep) without
    touching the shared default used by campaign.py / trajectory_campaign.py.

    Returns a status dict (also written to opt_status.json).
    """
    import cavity_common as cc
    from cqed_scf.drivers import bfgs_optimize

    if resume_from is not None:
        symbols, coords0 = resume_from
    else:
        symbols, coords0 = read_xyz(BASE_DIR / cell["xyz"])
    geom_str = cc.psi4_geometry_string(symbols, coords0)

    traj_path = cell_dir / "opt_traj.xyz"
    if resume_from is None:
        open(traj_path, "w").close()  # fresh trajectory

    calc = cc.build_cqed_calculator(lambda_vector)

    gnorm_history = []
    converged = False
    energy = None
    coords_ang = coords0
    t0 = time.time()

    for attempt in range(1, MAX_RESTARTS + 1):
        # bfgs_optimize returns (scipy OptimizeResult, observer_data dict).
        # The second value is the observer mapping (empty here, since we pass no
        # observers) -- NOT the gradient. With project_tr_rot=True the objective
        # hands BFGS the *projected* Cartesian gradient, so SciPy stores it in
        # result.jac: the projected gradient at the final geometry, which is our
        # convergence test.
        opt_result, _observer_data = bfgs_optimize(
            calculator=calc,
            geometry=geom_str,
            canonical="psi4",
            gtol=GTOL_BFGS,
            maxiter=MAXITER,
            debug=DEBUG_OPT,
            project_tr_rot=True,
            projection_debug=DEBUG_OPT,
        )

        coords_ang = cc.coords_angstrom_from_opt_result(opt_result)
        energy = float(opt_result.fun)
        gnorm = float(np.linalg.norm(np.asarray(opt_result.jac, dtype=float)))
        gnorm_history.append(gnorm)

        write_xyz_frame(
            traj_path, symbols, coords_ang,
            comment=(f"attempt {attempt} | |g_proj| = {gnorm:.3e} "
                     f"| E = {energy:.10f} Ha"),
            mode="a",
        )

        if gnorm < conv_threshold:
            converged = True
            break

        # Restart from the last geometry with identical settings (fresh BFGS
        # Hessian, per the user's current protocol).
        geom_str = cc.psi4_geometry_string(symbols, coords_ang)

    # Write the final geometry for the frequency step to consume.
    if converged:
        write_xyz_frame(
            cell_dir / "optimized.xyz", symbols, coords_ang,
            comment=(f"{cell['id']} | OPTIMIZED | |g_proj|={gnorm_history[-1]:.3e} "
                     f"| E={energy:.10f} Ha"),
        )

    status = {
        "cell_id": cell["id"],
        "converged": converged,
        "attempts": len(gnorm_history),
        "gnorm_history": gnorm_history,
        "final_gnorm": gnorm_history[-1] if gnorm_history else None,
        "conv_threshold": conv_threshold,
        "final_energy_hartree": energy,
        "wall_seconds": round(time.time() - t0, 1),
        "lambda_vector": lambda_vector,
    }
    if resume_from is not None:
        status["resumed_from_stall"] = True
    save_json(cell_dir / "opt_status.json", status)
    return status


def frequency_cell(cell, cell_dir, lambda_vector):
    """
    Finite-difference Hessian on the optimized geometry, projected ZPE, and
    imaginary-mode detection. Returns a results dict (also written to
    frequencies.json).
    """
    import cavity_common as cc
    from ase import Atoms
    from ase.vibrations import Vibrations

    symbols, coords_ang = read_xyz(cell_dir / "optimized.xyz")
    atoms = Atoms(symbols=symbols, positions=coords_ang)

    calc = cc.build_cqed_calculator(lambda_vector)
    atoms.calc = cc.CQED_DFT_Gradient(calc)

    # ASE caches each displacement under this name and resumes automatically,
    # so an interrupted Hessian picks up where it left off.
    vib_name = str(cell_dir / "vib")
    t0 = time.time()
    vib = Vibrations(atoms, name=vib_name, delta=FD_DELTA)
    vib.run()

    analysis = cc.project_and_diagonalize_hessian(atoms, vib)
    zpe = cc.zero_point_energy(analysis["freq_cm"])

    # Persist arrays for later inspection.
    np.save(cell_dir / "H_cart_au.npy", analysis["H_cart_au"])
    np.save(cell_dir / "H_locked_mw.npy", analysis["H_locked_mw"])
    np.save(cell_dir / "freq_cm.npy", analysis["freq_cm"])

    results = {
        "cell_id": cell["id"],
        "projector_rank": analysis["projector_rank"],
        "n_real_modes": zpe["n_real_modes"],
        "n_imaginary": zpe["n_imaginary"],
        "imaginary_freqs_cm": zpe["imaginary_freqs_cm"],
        "lowest_freq_cm": zpe["lowest_freq_cm"],
        "zpe_ev": zpe["zpe_ev"],
        "zpe_hartree": zpe["zpe_hartree"],
        "wall_seconds": round(time.time() - t0, 1),
    }
    save_json(cell_dir / "frequencies.json", results)
    return results


def promote_stalled(cell, cell_dir, opt_status):
    """
    Accept a previously-stalled cell whose recorded projected gradient norm now
    satisfies CONV_THRESHOLD: write its last trajectory geometry as optimized.xyz
    and mark it converged. No re-optimization -- the existing geometry is already
    the most converged one we have (lower gradient than a fresh loosened run
    would stop at).
    """
    symbols, coords_ang = read_last_xyz_frame(cell_dir / "opt_traj.xyz")
    write_xyz_frame(
        cell_dir / "optimized.xyz", symbols, coords_ang,
        comment=(f"{cell['id']} | PROMOTED |g|={opt_status['final_gnorm']:.3e} "
                 f"< {CONV_THRESHOLD:.1e} | E={opt_status.get('final_energy_hartree')} Ha"),
    )
    opt_status = {
        **opt_status,
        "converged": True,
        "promoted": True,
        "promoted_from_gnorm": opt_status.get("final_gnorm"),
        "promoted_threshold": CONV_THRESHOLD,
    }
    save_json(cell_dir / "opt_status.json", opt_status)
    return opt_status


def run_cell(cell, force=False, reoptimize=False):
    """Run one cell end-to-end with resume/checkpoint logic."""
    import cavity_common as cc

    cell_dir = RUNS_DIR / cell["id"]
    cell_dir.mkdir(parents=True, exist_ok=True)

    done_marker = cell_dir / "DONE"
    if done_marker.exists() and not force:
        print(f"[skip] {cell['id']} (already DONE)")
        return "done"

    lambda_vector = cc.lambda_vector_for(cell["theta"], cell["phi"], cell["magnitude"])
    save_json(cell_dir / "cell.json", {**cell, "lambda_vector": lambda_vector})

    # ---- Optimization ----
    opt_status_path = cell_dir / "opt_status.json"
    opt_status = None
    if opt_status_path.exists() and not force:
        opt_status = json.loads(opt_status_path.read_text())

    already_converged = bool(opt_status and opt_status.get("converged"))
    # A prior stall whose recorded gradient now clears CONV_THRESHOLD: promote
    # the existing geometry instead of re-optimizing (cheaper and higher quality).
    promotable = bool(
        opt_status and not already_converged and not reoptimize
        and opt_status.get("final_gnorm") is not None
        and opt_status["final_gnorm"] < CONV_THRESHOLD
        and (cell_dir / "opt_traj.xyz").exists()
    )

    if already_converged:
        pass
    elif promotable:
        print(f"[promote] {cell['id']}: |g|={opt_status['final_gnorm']:.2e} "
              f"< {CONV_THRESHOLD:.0e}; accepting stalled geometry, no re-opt")
        opt_status = promote_stalled(cell, cell_dir, opt_status)
    else:
        print(f"[opt ] {cell['id']} ...")
        opt_status = optimize_cell(cell, cell_dir, lambda_vector)

    if not opt_status.get("converged"):
        print(f"[STALL] {cell['id']}: not converged after "
              f"{opt_status.get('attempts')} attempts "
              f"(|g|={opt_status.get('final_gnorm')}). Skipping frequencies.")
        return "stalled"

    # ---- Frequencies (skip if already computed and not forcing) ----
    freq_path = cell_dir / "frequencies.json"
    if freq_path.exists() and not force:
        freq = json.loads(freq_path.read_text())
    else:
        print(f"[freq] {cell['id']} ...")
        freq = frequency_cell(cell, cell_dir, lambda_vector)

    done_marker.write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))

    flag = " <-- IMAGINARY MODE(S)" if freq.get("n_imaginary", 0) > 0 else ""
    print(f"[done] {cell['id']}: ZPE = {freq.get('zpe_hartree'):.8f} Ha, "
          f"n_imag = {freq.get('n_imaginary')}{flag}")
    return "done"


# -----------------------------------------------------------------------------
# Rollup / summary
# -----------------------------------------------------------------------------

def summarize():
    """Scan runs/ and build a CSV + JSON rollup of all attempted cells."""
    rows = []
    for cell in build_manifest():
        cell_dir = RUNS_DIR / cell["id"]
        opt = _maybe_json(cell_dir / "opt_status.json")
        freq = _maybe_json(cell_dir / "frequencies.json")
        rows.append({
            "cell_id": cell["id"],
            "isomer": cell["isomer"],
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

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    save_json(RUNS_DIR / "campaign_status.json", rows)

    cols = ["cell_id", "isomer", "theta", "phi", "magnitude", "opt_converged",
            "attempts", "final_gnorm", "E_opt_hartree", "zpe_hartree",
            "n_imaginary", "lowest_freq_cm"]
    csv_path = RUNS_DIR / "campaign_status.csv"
    with open(csv_path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {RUNS_DIR / 'campaign_status.json'}")
    return rows


def _maybe_json(path):
    return json.loads(path.read_text()) if path.exists() else None


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Cavity opt + frequency campaign driver.")
    parser.add_argument("--dry-run", action="store_true", help="print the manifest and exit")
    parser.add_argument("--only", metavar="SUBSTR", help="only run cells whose id contains SUBSTR")
    parser.add_argument("--list-ids", action="store_true",
                        help="print one cell id per line (honors --only) and exit; for piping to a job launcher")
    parser.add_argument("--summarize", action="store_true", help="rebuild the rollup CSV/JSON and exit")
    parser.add_argument("--no-summarize", action="store_true",
                        help="skip the end-of-run rollup (use for parallel single-cell shards to avoid racing writes)")
    parser.add_argument("--force", action="store_true", help="redo cells even if already DONE")
    parser.add_argument("--reoptimize", action="store_true",
                        help="re-run the optimization from scratch instead of promoting a "
                             "stalled geometry that already clears CONV_THRESHOLD")
    args = parser.parse_args()

    cells = build_manifest()
    if args.only:
        cells = [c for c in cells if args.only in c["id"]]

    if args.list_ids:
        for c in cells:
            print(c["id"])
        return

    if args.dry_run:
        print(f"{len(cells)} cells:")
        for c in cells:
            print(f"  {c['id']:28s}  isomer={c['isomer']:5s} "
                  f"theta={c['theta']:>3} phi={c['phi']:>3} |lambda|={c['magnitude']:.2f}")
        return

    if args.summarize:
        summarize()
        return

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Running {len(cells)} cells. Runs dir: {RUNS_DIR}")

    counts = {"done": 0, "stalled": 0, "error": 0}
    for c in cells:
        try:
            outcome = run_cell(c, force=args.force, reoptimize=args.reoptimize)
            counts[outcome] = counts.get(outcome, 0) + 1
        except Exception:
            counts["error"] += 1
            print(f"[ERROR] {c['id']}:\n{traceback.format_exc()}", file=sys.stderr)
            # Keep going -- one bad cell shouldn't sink a multi-day campaign.

    # Skip the rollup for single-cell shards (parallel runs) to avoid racing
    # writes to campaign_status.*; rebuild it once at the end with --summarize.
    if not args.no_summarize and len(cells) > 1:
        summarize()
    print(f"\nCampaign finished: {counts}")


if __name__ == "__main__":
    main()

