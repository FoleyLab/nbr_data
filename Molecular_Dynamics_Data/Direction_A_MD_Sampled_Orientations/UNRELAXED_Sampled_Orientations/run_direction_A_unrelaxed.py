#!/usr/bin/env python3
"""Run unrelaxed Direction A CQED-DFT single-point calculations.

This script combines the 25 sampled cavity orientations in
trajectory_sampling_direction_A.dat with the meta, para, and ortho unrelaxed
XYZ structures, producing 75 single-point energy calculations.
"""

from __future__ import annotations

import argparse
import csv
import os
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SAMPLING_FILE = ROOT / "trajectory_sampling_direction_A.dat"
ISOMER_XYZ = {
    "meta": ROOT / "unrelaxed_meta.xyz",
    "para": ROOT / "unrelaxed_para.xyz",
    "ortho": ROOT / "unrelaxed_ortho.xyz",
}

PSI4_OPTIONS = {
    "basis": "6-311G*",
    "scf_type": "df",
    "e_convergence": 1e-12,
    "d_convergence": 1e-12,
    "dft_radial_points": 99,
    "dft_spherical_points": 590,
    "dft_pruning_scheme": "none",
}

FIELDNAMES = [
    "job_id",
    "isomer",
    "sample_index",
    "theta_deg",
    "phi_deg",
    "lambda_mag",
    "lambda_x",
    "lambda_y",
    "lambda_z",
    "omega",
    "charge",
    "multiplicity",
    "functional",
    "basis",
    "energy_hartree",
    "status",
    "log_file",
    "psi4_output_file",
]


def parse_sampling_file(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip()),
            delimiter="\t",
        )
        if reader.fieldnames and reader.fieldnames[0].startswith("#"):
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")

        for index, row in enumerate(reader, start=1):
            rows.append(
                {
                    "sample_index": f"{index:03d}",
                    "theta_deg": row["Theta"],
                    "phi_deg": row["Phi"],
                    "lambda_mag": row["lambda_mag"],
                    "omega": row["omega"],
                    "charge": row["charge"],
                    "multiplicity": row["multiplicity"],
                    "functional": row["functional"],
                    "basis": row["basis"],
                }
            )
    return rows


def xyz_to_psi4_geometry(path: Path, charge: int, multiplicity: int) -> str:
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path} does not look like a valid XYZ file.")

    try:
        atom_count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"First line of {path} must be the atom count.") from exc

    atom_lines = [line.rstrip() for line in lines[2 : 2 + atom_count]]
    if len(atom_lines) != atom_count:
        raise ValueError(f"{path} ended before reading {atom_count} atoms.")

    return "\n".join(
        [
            f"{charge} {multiplicity}",
            *atom_lines,
            "units angstrom",
            "no_reorient",
            "no_com",
            "symmetry c1",
            "",
        ]
    )


def build_jobs(isomers: list[str]) -> list[dict[str, object]]:
    samples = parse_sampling_file(SAMPLING_FILE)
    jobs: list[dict[str, object]] = []

    for isomer in isomers:
        for sample in samples:
            charge = int(sample["charge"])
            multiplicity = int(sample["multiplicity"])
            theta = int(float(sample["theta_deg"]))
            phi = int(float(sample["phi_deg"]))
            sample_index = sample["sample_index"]
            job_id = f"{isomer}_dirA_r{sample_index}_th{theta}_ph{phi}"

            jobs.append(
                {
                    "job_id": job_id,
                    "isomer": isomer,
                    "sample_index": sample_index,
                    "theta_deg": theta,
                    "phi_deg": phi,
                    "lambda_mag": float(sample["lambda_mag"]),
                    "omega": float(sample["omega"]),
                    "charge": charge,
                    "multiplicity": multiplicity,
                    "functional": sample["functional"],
                    "basis": sample["basis"],
                    "geometry": xyz_to_psi4_geometry(
                        ISOMER_XYZ[isomer],
                        charge=charge,
                        multiplicity=multiplicity,
                    ),
                }
            )

    return jobs


def configure_worker(memory: str, threads: int, scratch_root: str) -> None:
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)

    import psi4

    psi4.core.be_quiet()
    psi4.set_memory(memory)
    psi4.set_num_threads(threads)
    psi4.set_options(PSI4_OPTIONS)

    scratch = Path(scratch_root) / f"psi4_direction_A_worker_{os.getpid()}"
    scratch.mkdir(parents=True, exist_ok=True)
    psi4.core.IOManager.shared_object().set_default_path(str(scratch))


def run_job(job: dict[str, object], output_dir: str) -> dict[str, object]:
    output_path = Path(output_dir)
    logs_dir = output_path / "logs"
    psi4_dir = output_path / "psi4_outputs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    psi4_dir.mkdir(parents=True, exist_ok=True)

    job_id = str(job["job_id"])
    log_file = logs_dir / f"{job_id}.log"
    psi4_output_file = psi4_dir / f"{job_id}.dat"

    row = {
        "job_id": job_id,
        "isomer": job["isomer"],
        "sample_index": job["sample_index"],
        "theta_deg": job["theta_deg"],
        "phi_deg": job["phi_deg"],
        "lambda_mag": job["lambda_mag"],
        "lambda_x": "",
        "lambda_y": "",
        "lambda_z": "",
        "omega": job["omega"],
        "charge": job["charge"],
        "multiplicity": job["multiplicity"],
        "functional": job["functional"],
        "basis": job["basis"],
        "energy_hartree": "",
        "status": "",
        "log_file": str(log_file),
        "psi4_output_file": str(psi4_output_file),
    }

    with log_file.open("w") as log_handle, redirect_stdout(log_handle), redirect_stderr(log_handle):
        try:
            import numpy as np
            import psi4
            from cqed_scf import CQEDCalculator
            from cqed_scf.utils import generate_field_vector_from_theta_and_phi

            psi4.core.set_output_file(str(psi4_output_file), False)
            psi4.core.clean()
            psi4.core.clean_options()
            psi4.set_options(PSI4_OPTIONS)

            theta = float(job["theta_deg"])
            phi = float(job["phi_deg"])
            lambda_mag = float(job["lambda_mag"])

            field_vector = generate_field_vector_from_theta_and_phi(theta, phi)
            lambda_direction = np.asarray(field_vector, dtype=float)
            lambda_direction /= np.linalg.norm(lambda_direction)
            lambda_vector = (lambda_mag * lambda_direction).tolist()

            row["lambda_x"], row["lambda_y"], row["lambda_z"] = lambda_vector

            print(f"Job: {job_id}")
            print(f"Isomer: {job['isomer']}")
            print(f"Sample index: {job['sample_index']}")
            print(f"theta={theta:g} phi={phi:g}")
            print(f"|lambda|={lambda_mag:g} lambda_vector={lambda_vector}")
            print(f"omega={job['omega']}")
            print(f"functional={job['functional']} basis={job['basis']}")

            calc = CQEDCalculator(
                lambda_vector=lambda_vector,
                psi4_options=PSI4_OPTIONS,
                omega=float(job["omega"]),
                density_fitting=True,
                charge=int(job["charge"]),
                multiplicity=int(job["multiplicity"]),
                functional=str(job["functional"]),
                debug=False,
            )
            energy = float(calc.energy(str(job["geometry"])))

            row["energy_hartree"] = f"{energy:.16f}"
            row["status"] = "ok"
            print(f"Energy: {energy:.16f} Hartree")
            psi4.core.clean()
        except Exception as exc:
            row["status"] = (
                f"failed: {type(exc).__name__}: {exc}; "
                f"traceback: {traceback.format_exc().replace(os.linesep, ' | ')}"
            )

    return row


def load_finished_job_ids(output_csv: Path) -> set[str]:
    if not output_csv.exists():
        return set()

    finished: set[str] = set()
    with output_csv.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("status") == "ok":
                finished.add(row["job_id"])
    return finished


def append_row(output_csv: Path, row: dict[str, object]) -> None:
    exists = output_csv.exists()
    with output_csv.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run 75 unrelaxed Direction A CQED-DFT single-point calculations."
    )
    parser.add_argument("--workers", type=int, default=9, help="Concurrent calculations.")
    parser.add_argument("--threads-per-worker", type=int, default=1)
    parser.add_argument("--memory-per-worker", default="2 GB")
    parser.add_argument("--output-dir", default=str(ROOT / "direction_A_unrelaxed_results"))
    parser.add_argument(
        "--scratch-root",
        default="/tmp",
        help="Root directory for per-worker Psi4 scratch folders.",
    )
    parser.add_argument(
        "--isomers",
        nargs="+",
        choices=sorted(ISOMER_XYZ),
        default=["meta", "para", "ortho"],
        help="Subset of isomers to run.",
    )
    parser.add_argument("--dry-run", action="store_true", help="List jobs without running Psi4.")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_csv = output_dir / "direction_A_unrelaxed_energies.csv"
    jobs = build_jobs(args.isomers)
    finished = load_finished_job_ids(output_csv)
    pending = [job for job in jobs if str(job["job_id"]) not in finished]

    print(f"Sampling file: {SAMPLING_FILE}")
    print(f"Output directory: {output_dir}")
    print(f"Total jobs selected: {len(jobs)}")
    print(f"Already completed: {len(finished)}")
    print(f"Pending jobs: {len(pending)}")
    print(f"Workers: {args.workers}")
    print(f"Threads per worker: {args.threads_per_worker}")
    print(f"Memory per worker: {args.memory_per_worker}")

    if args.dry_run:
        for job in pending:
            print(
                f"{job['job_id']}: isomer={job['isomer']} "
                f"theta={job['theta_deg']} phi={job['phi_deg']} "
                f"lambda={job['lambda_mag']} omega={job['omega']}"
            )
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    if args.workers == 1:
        configure_worker(args.memory_per_worker, args.threads_per_worker, args.scratch_root)
        for job in pending:
            row = run_job(job, str(output_dir))
            append_row(output_csv, row)
            print(f"{row['job_id']}: {row['status']} energy={row['energy_hartree']}")
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=configure_worker,
            initargs=(args.memory_per_worker, args.threads_per_worker, args.scratch_root),
        ) as executor:
            futures = [executor.submit(run_job, job, str(output_dir)) for job in pending]
            for future in as_completed(futures):
                row = future.result()
                append_row(output_csv, row)
                print(f"{row['job_id']}: {row['status']} energy={row['energy_hartree']}")

    print(f"Results written to {output_csv}")


if __name__ == "__main__":
    main()
