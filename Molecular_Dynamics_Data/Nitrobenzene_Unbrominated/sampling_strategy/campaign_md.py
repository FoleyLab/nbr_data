#!/usr/bin/env python3
"""
campaign_md.py -- trajectory-aware driver for the MD orientation sweep.
=======================================================================

Sits alongside the grid `campaign.py`. Instead of iterating (isomer, direction,
magnitude) cells, it iterates cells defined by manifest.json:

    cell = (intermediate, orientation.id)          intermediate in {ortho,meta,para}

For each cell it reproduces the SAME per-cell workflow as the grid campaign:

    derive lambda once  ->  constrained optimize (restart loop gated on the
    PROJECTED gradient norm < CONV_THRESHOLD)  ->  gate  ->  ASE finite-diff
    Hessian on the FULL gradient  ->  rigid-body-projected ZPE  ->  archive in
    runs/<cell_id>/ with a DONE marker.

The only genuinely new step is deriving lambda: rather than
lambda_vector_for(theta, phi, mag) in the lab frame, we transfer the manifest's
body-frame vector (a,b,c) onto THIS intermediate's own geometry via
orientation_transfer.lambda_for_intermediate (which asserts Br-on-minus-z and
frame orthonormality). Everything downstream is unchanged, so resumability and
`--list-ids | parallel ... --only {} --no-summarize` carry over verbatim.

Assumed interfaces (from cavity_common.py, per the project primer)
------------------------------------------------------------------
    build_cqed_calculator(geometry_string, lambda_vec, **cfg) -> calculator
    psi4_geometry_string(symbols, coords, charge, multiplicity) -> str
    CQED_DFT_Gradient(calculator, ...)                 # ASE calculator
    bfgs_optimize(atoms, project_tr_rot=True, ...) -> scipy-like result (.x, .jac)
    project_and_diagonalize_hessian(H, masses, coords) -> (freqs_cm, modes)
    zero_point_energy(freqs_cm) -> zpe_hartree
    rigid_body_projector_mw(masses, coords)            # used inside the above

Intermediate reference geometries
----------------------------------
    intermediates/<name>.xyz   (charge +1, singlet bromo-nitrobenzene cation)
Atom ordering must match the reactant frame convention (ipso C nearest N;
first three carbons define the ring normal). If a given intermediate lists atoms
differently, set an explicit triplet in RING_TRIPLET_OVERRIDE below; otherwise
the Br-on-minus-z assert will stop the run before any CPU is wasted.
"""

import argparse, json, os, sys, time
import numpy as np

from orientation_transfer import lambda_for_intermediate

# ------------------------------------------------------------------- configuration
CONV_THRESHOLD = 5e-4               # projected-gradient norm gate (same as grid)
MAX_RESTARTS   = 8
CHARGE, MULT   = 1, 1               # bromo-nitrobenzene cation, singlet
RUNS_DIR       = "runs"
INTER_DIR      = "intermediates"
CQED_CFG = dict(method="wb97x", basis="6-311G*", density_fit=True, omega=0.06615)

# Per-intermediate ring-triplet override (None => first three carbons in file
# order). Fill in only if an intermediate's atom ordering differs; e.g.
#   RING_TRIPLET_OVERRIDE = {"meta": [0, 1, 3]}
RING_TRIPLET_OVERRIDE = {}


# --------------------------------------------------------------------------- io util
def load_manifest(path):
    with open(path) as fh:
        return json.load(fh)


def read_xyz(path):
    lines = open(path).read().splitlines()
    nat = int(lines[0].strip())
    syms, xyz = [], []
    for j in range(nat):
        p = lines[2 + j].split()
        syms.append(p[0]); xyz.append([float(x) for x in p[1:4]])
    return syms, np.array(xyz)


def write_xyz(path, symbols, coords, comment=""):
    with open(path, "w") as fh:
        fh.write(f"{len(symbols)}\n{comment}\n")
        for s, r in zip(symbols, coords):
            fh.write(f"{s:2s} {r[0]:18.10f} {r[1]:18.10f} {r[2]:18.10f}\n")


def cell_ids(manifest):
    for o in manifest["orientations"]:
        for inter in manifest["meta"]["intermediates"]:
            yield f"{inter}__{o['id']}", inter, o


# ----------------------------------------------------------------------- per-cell run
def run_cell(cell_id, inter, orient, mag, force=False):
    cell_dir = os.path.join(RUNS_DIR, cell_id)
    done = os.path.join(cell_dir, "DONE")
    if os.path.exists(done) and not force:
        return "skip (DONE)"
    os.makedirs(cell_dir, exist_ok=True)

    # ---- heavy imports deferred so --list-ids etc. work without Psi4 present
    from cavity_common import (build_cqed_calculator, psi4_geometry_string,
                               CQED_DFT_Gradient, bfgs_optimize,
                               project_and_diagonalize_hessian, zero_point_energy)
    from ase import Atoms
    from ase.vibrations import Vibrations

    symbols, coords0 = read_xyz(os.path.join(INTER_DIR, f"{inter}.xyz"))
    triplet = RING_TRIPLET_OVERRIDE.get(inter)

    # ---- derive lambda ONCE for this cell (asserts Br-on-minus-z, orthonormality)
    lam = lambda_for_intermediate(orient["abc"], symbols, coords0, mag,
                                  ring_triplet=triplet)
    np.save(os.path.join(cell_dir, "lambda.npy"), lam)

    # ---- constrained optimization with restart loop gated on projected gnorm
    coords = coords0.copy()
    converged, gnorm, restarts = False, np.inf, 0
    while restarts < MAX_RESTARTS:
        geom = psi4_geometry_string(symbols, coords, CHARGE, MULT)
        calc = build_cqed_calculator(geom, lam, **CQED_CFG)
        atoms = Atoms(symbols=symbols, positions=coords)
        atoms.calc = CQED_DFT_Gradient(calc)
        result = bfgs_optimize(atoms, project_tr_rot=True)
        coords = np.asarray(result.x).reshape(-1, 3)
        gnorm = float(np.linalg.norm(result.jac))     # PROJECTED gradient lives in .jac
        restarts += 1
        if gnorm < CONV_THRESHOLD:
            converged = True
            break

    # ---- gate: promotion of stalled-but-already-converged cells
    promoted = (not converged) and (gnorm < CONV_THRESHOLD)
    if not (converged or promoted):
        json.dump({"converged": False, "final_gnorm": gnorm, "restarts": restarts,
                   "reason": "gate_failed"},
                  open(os.path.join(cell_dir, "opt_status.json"), "w"), indent=2)
        return f"FAIL gate (gnorm={gnorm:.2e})"

    final_energy = float(atoms.get_potential_energy())
    write_xyz(os.path.join(cell_dir, "optimized.xyz"), symbols, coords,
              comment=f"{cell_id} E={final_energy:.10f} gnorm={gnorm:.2e}")
    json.dump({"final_energy_hartree": final_energy, "converged": bool(converged),
               "promoted": bool(promoted), "final_gnorm": gnorm,
               "restarts": restarts, "lambda": lam.tolist(),
               "abc": orient["abc"], "intermediate": inter},
              open(os.path.join(cell_dir, "opt_status.json"), "w"), indent=2)

    # ---- ASE finite-difference Hessian on the FULL gradient, then projected ZPE
    vib = Vibrations(atoms, name=os.path.join(cell_dir, "vib"))
    vib.run()
    H = vib.get_vibrations().get_hessian_2d()          # full (3N,3N) Hessian
    np.save(os.path.join(cell_dir, "hessian.npy"), H)
    masses = atoms.get_masses()
    freqs_cm, _ = project_and_diagonalize_hessian(H, masses, coords)
    zpe = zero_point_energy(freqs_cm)
    n_imag = int(np.sum(np.asarray(freqs_cm) < 0))
    lowest = float(np.min(np.abs(freqs_cm)))
    json.dump({"zpe_hartree": float(zpe), "n_imaginary": n_imag,
               "lowest_freq_cm": lowest,
               "final_energy_hartree": final_energy,
               "e_plus_zpe_hartree": final_energy + float(zpe)},
              open(os.path.join(cell_dir, "frequencies.json"), "w"), indent=2)

    open(done, "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    flag = "promoted" if promoted else "converged"
    return f"OK ({flag}) E+ZPE={final_energy + zpe:.6f} n_imag={n_imag}"


# ----------------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="manifest.json")
    ap.add_argument("--list-ids", action="store_true")
    ap.add_argument("--only", help="run a single cell id")
    ap.add_argument("--force", action="store_true", help="ignore DONE markers")
    ap.add_argument("--no-summarize", action="store_true")
    args = ap.parse_args()

    manifest = load_manifest(args.manifest)
    mag = manifest["meta"]["lambda_magnitude"]
    cells = list(cell_ids(manifest))

    if args.list_ids:
        for cid, _, _ in cells:
            print(cid)
        return

    targets = [c for c in cells if (args.only is None or c[0] == args.only)]
    if args.only and not targets:
        sys.exit(f"unknown cell id: {args.only}")

    results = {}
    for cid, inter, orient in targets:
        t0 = time.time()
        try:
            msg = run_cell(cid, inter, orient, mag, force=args.force)
        except AssertionError as e:              # frame mismatch -> loud stop
            msg = f"ASSERT {e}"
        except Exception as e:
            msg = f"ERROR {type(e).__name__}: {e}"
        results[cid] = msg
        print(f"[{cid}] {msg}  ({time.time() - t0:.0f}s)")

    if not args.no_summarize:
        ok = sum(v.startswith("OK") for v in results.values())
        print(f"\n{ok}/{len(results)} cells OK")


if __name__ == "__main__":
    main()
