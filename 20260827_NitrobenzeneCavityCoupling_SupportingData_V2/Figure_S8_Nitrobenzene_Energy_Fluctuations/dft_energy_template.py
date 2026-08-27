"""
Compute cavity-free single-point energies for each frame in an MD trajectory.

The trajectory already stores the cavity energy in each XYZ comment line.  This
driver reads those values, recomputes the same geometry with zero cavity field,
and writes both energies plus their difference to CSV.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import psi4

from cqed_scf import CQEDCalculator, CQEDConfig


HARTREE_TO_KCAL_MOL = 627.5094740631
THIS_DIR = Path(__file__).resolve().parent
DEFAULT_TRAJECTORY = THIS_DIR / "nitrobenzene_direction_A_wb97x_d_4000_ts.xyz"
DEFAULT_OUTPUT = THIS_DIR / "nitrobenzene_direction_A_wb97x_d_cavity_free_energies.csv"

ENERGY_RE = re.compile(r"\bE=([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)")

PSI4_OPTIONS = {
    "basis": "6-311g*",
    "scf_type": "df",
    "e_convergence": 1e-10,
    "d_convergence": 1e-9,
    "dft_radial_points": 99,
    "dft_spherical_points": 590,
    "dft_pruning_scheme": "none",
}

NO_CAVITY_FIELD_VECTOR = np.array([0.0, 0.0, 0.0])
NO_CAVITY_OMEGA = 0.0

_CALC_NO_CAVITY = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute zero-field wb97x-d energies for frames in a CQED-MD XYZ trajectory."
    )
    parser.add_argument("--trajectory", type=Path, default=DEFAULT_TRAJECTORY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--nproc", type=int, default=1, help="Number of frame calculations to run at once.")
    parser.add_argument("--threads-per-worker", type=int, default=1, help="Psi4 threads used inside each worker.")
    parser.add_argument("--memory", default="4 GB", help="Psi4 memory per worker.")
    parser.add_argument("--start", type=int, default=0, help="First frame index to process.")
    parser.add_argument("--stop", type=int, default=None, help="Stop before this frame index.")
    parser.add_argument("--stride", type=int, default=1, help="Process every Nth frame starting from --start.")
    parser.add_argument("--resume", action="store_true", help="Skip frame indices already present in the CSV.")
    args = parser.parse_args()
    if args.stride < 1:
        parser.error("--stride must be at least 1")
    return args


def iter_xyz_frames(filename):
    with open(filename) as handle:
        frame_index = 0
        while True:
            natom_line = handle.readline()
            if not natom_line:
                break

            natom_line = natom_line.strip()
            if not natom_line:
                continue

            natom = int(natom_line)
            comment = handle.readline().strip()
            atoms = [handle.readline().rstrip() for _ in range(natom)]

            match = ENERGY_RE.search(comment)
            if match is None:
                raise ValueError(f"Could not parse cavity energy from frame {frame_index}: {comment}")

            yield {
                "frame": frame_index,
                "comment": comment,
                "cavity_energy_hartree": float(match.group(1)),
                "geometry": build_psi4_geometry(atoms),
            }
            frame_index += 1


def build_psi4_geometry(atom_lines):
    return "\n".join(
        [
            "0 1",
            *atom_lines,
            "units angstrom",
            "no_reorient",
            "symmetry c1",
        ]
    )


def completed_frames(csv_file):
    if not csv_file.exists():
        return set()

    with open(csv_file, newline="") as handle:
        reader = csv.DictReader(handle)
        return {int(row["frame"]) for row in reader if row.get("frame")}


def initialize_worker(memory, threads_per_worker):
    global _CALC_NO_CAVITY

    psi4.core.be_quiet()
    psi4.set_memory(memory)
    psi4.set_num_threads(threads_per_worker)
    psi4.core.set_output_file(f"psi4_no_cavity_worker_{os.getpid()}.out", False)

    config_no_cavity = CQEDConfig(
        lambda_vector=NO_CAVITY_FIELD_VECTOR,
        omega=NO_CAVITY_OMEGA,
        psi4_options=PSI4_OPTIONS,
        reference="rks",
        functional="wb97x-d",
        density_fitting=True,
        charge=0,
        multiplicity=1,
        dispersion_policy="post_scf",
        debug=False,
        quiet=True,
    )
    _CALC_NO_CAVITY = CQEDCalculator(config=config_no_cavity)


def compute_no_cavity_energy(frame):
    psi4.core.clean()
    psi4.core.clean_options()

    no_cavity_energy = _CALC_NO_CAVITY.energy(frame["geometry"])
    energy_difference = frame["cavity_energy_hartree"] - no_cavity_energy

    return {
        "frame": frame["frame"],
        "cavity_energy_hartree": frame["cavity_energy_hartree"],
        "no_cavity_energy_hartree": no_cavity_energy,
        "delta_cavity_minus_no_cavity_hartree": energy_difference,
        "delta_cavity_minus_no_cavity_kcal_mol": energy_difference * HARTREE_TO_KCAL_MOL,
        "xyz_comment": frame["comment"],
    }


def write_rows_in_order(output_file, rows, append):
    fieldnames = [
        "frame",
        "cavity_energy_hartree",
        "no_cavity_energy_hartree",
        "delta_cavity_minus_no_cavity_hartree",
        "delta_cavity_minus_no_cavity_kcal_mol",
        "xyz_comment",
    ]

    with open(output_file, "a" if append else "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if not append:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)
        handle.flush()


def main():
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    skip_frames = completed_frames(args.output) if args.resume else set()
    frames = [
        frame
        for frame in iter_xyz_frames(args.trajectory)
        if frame["frame"] >= args.start
        and (args.stop is None or frame["frame"] < args.stop)
        and (frame["frame"] - args.start) % args.stride == 0
        and frame["frame"] not in skip_frames
    ]

    if not frames:
        print("No frames to process.")
        return

    print(f"Trajectory : {args.trajectory}")
    print(f"Output     : {args.output}")
    print(f"Frames     : {frames[0]['frame']} through {frames[-1]['frame']} ({len(frames)} total)")
    print(f"Stride     : every {args.stride} frame(s)")
    print(f"Workers    : {args.nproc}")
    print(f"Psi4/thread: {args.threads_per_worker} per worker")
    print(f"No cavity  : lambda={NO_CAVITY_FIELD_VECTOR}, omega={NO_CAVITY_OMEGA}")

    append = args.resume and args.output.exists()
    pending = {}
    frame_order = [frame["frame"] for frame in frames]
    next_frame_order_index = 0

    with ProcessPoolExecutor(
        max_workers=args.nproc,
        initializer=initialize_worker,
        initargs=(args.memory, args.threads_per_worker),
    ) as executor:
        futures = {executor.submit(compute_no_cavity_energy, frame): frame["frame"] for frame in frames}

        for future in as_completed(futures):
            frame_index = futures[future]
            row = future.result()
            pending[frame_index] = row

            ready_rows = []
            while (
                next_frame_order_index < len(frame_order)
                and frame_order[next_frame_order_index] in pending
            ):
                next_frame_to_write = frame_order[next_frame_order_index]
                ready_rows.append(pending.pop(next_frame_to_write))
                next_frame_order_index += 1

            if ready_rows:
                write_rows_in_order(args.output, ready_rows, append=append)
                append = True
                last = ready_rows[-1]
                print(
                    f"Wrote frame {last['frame']:5d}: "
                    f"dE = {last['delta_cavity_minus_no_cavity_kcal_mol']: .8f} kcal/mol"
                )


if __name__ == "__main__":
    main()
