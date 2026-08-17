#!/usr/bin/env python3
"""Analyze full vs projected CQED-DFT gradients on theta/phi orientation rings.

This script descends into folders named like ``ortho_th72_ph180``, reads
``optimized.xyz`` and ``cell.json``, computes full/projected gradients, and writes
per-isomer plus pairwise diagnostics for testing torque-zero hypotheses.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


DEFAULT_ROOT = Path(__file__).resolve().parent
DEFAULT_ISOMERS = ("ortho", "meta")
DEFAULT_PSI4_OPTIONS = {
    "basis": "6-311G*",
    "scf_type": "df",
    "e_convergence": 1e-10,
    "d_convergence": 1e-8,
}
AU_TO_KCAL = 627.509


@dataclass(frozen=True)
class OrientationCase:
    isomer: str
    theta: float
    phi: float
    folder: Path
    geometry_path: Path
    cell_path: Path
    lambda_vector: np.ndarray


def parse_numeric_list(spec: str) -> list[float]:
    """Parse comma values or inclusive start:stop:step syntax."""
    if ":" in spec:
        parts = [float(part) for part in spec.split(":")]
        if len(parts) != 3:
            raise ValueError(f"Expected start:stop:step, got {spec!r}")
        start, stop, step = parts
        if step == 0:
            raise ValueError("Step cannot be zero")
        values = []
        value = start
        if step > 0:
            while value <= stop + abs(step) * 1e-9:
                values.append(round(value, 10))
                value += step
        else:
            while value >= stop - abs(step) * 1e-9:
                values.append(round(value, 10))
                value += step
        return values

    return [float(part.strip()) for part in spec.split(",") if part.strip()]


def label_number(value: float) -> str:
    if math.isclose(value, round(value), abs_tol=1e-8):
        return str(int(round(value)))
    return f"{value:g}".replace(".", "p")


def folder_for(root: Path, isomer: str, theta: float, phi: float) -> Path:
    return root / f"{isomer}_th{label_number(theta)}_ph{label_number(phi)}"


def read_lambda_vector(cell_path: Path) -> np.ndarray:
    with cell_path.open() as handle:
        cell = json.load(handle)
    try:
        return np.array(cell["lambda_vector"], dtype=float)
    except KeyError as exc:
        raise KeyError(f"{cell_path} does not contain cell['lambda_vector']") from exc


def xyz_to_psi4_geometry(xyz_path: Path, charge: int, multiplicity: int) -> str:
    lines = xyz_path.read_text().splitlines()
    if not lines:
        raise ValueError(f"{xyz_path} is empty")

    first = lines[0].strip()
    atom_lines = lines[2:] if first.isdigit() else lines
    atom_lines = [line.rstrip() for line in atom_lines if line.strip()]

    return "\n".join(
        [
            f"{charge} {multiplicity}",
            *atom_lines,
            "no_reorient",
            "no_com",
            "symmetry c1",
            "",
        ]
    )


def discover_cases(
    root: Path,
    isomers: Iterable[str],
    thetas: Iterable[float],
    phis: Iterable[float],
) -> list[OrientationCase]:
    cases = []
    for theta in thetas:
        for phi in phis:
            for isomer in isomers:
                folder = folder_for(root, isomer, theta, phi)
                geometry_path = folder / "optimized.xyz"
                cell_path = folder / "cell.json"
                if not folder.is_dir():
                    raise FileNotFoundError(f"Missing folder: {folder}")
                if not geometry_path.is_file():
                    raise FileNotFoundError(f"Missing geometry: {geometry_path}")
                if not cell_path.is_file():
                    raise FileNotFoundError(f"Missing cell metadata: {cell_path}")

                cases.append(
                    OrientationCase(
                        isomer=isomer,
                        theta=theta,
                        phi=phi,
                        folder=folder,
                        geometry_path=geometry_path,
                        cell_path=cell_path,
                        lambda_vector=read_lambda_vector(cell_path),
                    )
                )
    return cases


def warn_if_lambda_mismatch(cases: list[OrientationCase], tolerance: float = 1e-10) -> None:
    by_orientation: dict[tuple[float, float], list[OrientationCase]] = {}
    for case in cases:
        by_orientation.setdefault((case.theta, case.phi), []).append(case)

    for (theta, phi), orientation_cases in by_orientation.items():
        reference = orientation_cases[0].lambda_vector
        for case in orientation_cases[1:]:
            delta = np.linalg.norm(case.lambda_vector - reference)
            if delta > tolerance:
                print(
                    f"WARNING: lambda mismatch at theta={theta:g}, phi={phi:g}; "
                    f"{orientation_cases[0].isomer} vs {case.isomer} delta={delta:.3e}",
                    file=sys.stderr,
                )


def as_numpy_gradient(gradient) -> np.ndarray:
    if hasattr(gradient, "to_array"):
        gradient = gradient.to_array()
    return np.asarray(gradient, dtype=float)


def make_calculator(lambda_vector: np.ndarray, args):
    from cqed_scf import CQEDCalculator, CQEDConfig

    psi4_options = {
        "basis": args.basis,
        "scf_type": args.scf_type,
        "e_convergence": args.e_convergence,
        "d_convergence": args.d_convergence,
    }

    config = CQEDConfig(
        lambda_vector=lambda_vector,
        omega=args.omega,
        psi4_options=psi4_options,
        reference=args.reference,
        functional=args.functional,
        density_fitting=args.density_fitting,
        charge=args.charge,
        multiplicity=args.multiplicity,
        dispersion_policy=args.dispersion_policy,
        debug=args.debug,
        quiet=args.quiet,
    )
    return CQEDCalculator(config=config)


def compute_case(case: OrientationCase, args) -> dict[str, object]:
    geometry = xyz_to_psi4_geometry(case.geometry_path, args.charge, args.multiplicity)
    calc = make_calculator(case.lambda_vector, args)

    energy_full, gradient_full, _ = calc.energy_and_gradient(geometry)
    energy_projected, gradient_projected, _ = calc.energy_and_projected_gradient(geometry)

    gradient_full = as_numpy_gradient(gradient_full)
    gradient_projected = as_numpy_gradient(gradient_projected)
    gradient_removed = gradient_full - gradient_projected

    norm_full = float(np.linalg.norm(gradient_full))
    norm_projected = float(np.linalg.norm(gradient_projected))
    norm_removed = float(np.linalg.norm(gradient_removed))

    return {
        "theta": case.theta,
        "phi": case.phi,
        "isomer": case.isomer,
        "energy_full_hartree": float(energy_full),
        "energy_projected_hartree": float(energy_projected),
        "delta_energy_projected_minus_full_hartree": float(energy_projected - energy_full),
        "norm_full": norm_full,
        "norm_projected": norm_projected,
        "norm_removed": norm_removed,
        "lambda_x": float(case.lambda_vector[0]),
        "lambda_y": float(case.lambda_vector[1]),
        "lambda_z": float(case.lambda_vector[2]),
        "geometry_path": str(case.geometry_path),
    }


def write_csv(rows: list[dict[str, object]], path: Path, fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_pairwise_rows(
    rows: list[dict[str, object]], isomer_a: str, isomer_b: str
) -> list[dict[str, object]]:
    by_key = {(row["theta"], row["phi"], row["isomer"]): row for row in rows}
    orientations = sorted({(row["theta"], row["phi"]) for row in rows})
    pair_rows = []

    for theta, phi in orientations:
        row_a = by_key.get((theta, phi, isomer_a))
        row_b = by_key.get((theta, phi, isomer_b))
        if row_a is None or row_b is None:
            continue

        e_a = float(row_a["energy_full_hartree"])
        e_b = float(row_b["energy_full_hartree"])
        removed_a = float(row_a["norm_removed"])
        removed_b = float(row_b["norm_removed"])

        pair_rows.append(
            {
                "theta": theta,
                "phi": phi,
                "isomer_a": isomer_a,
                "isomer_b": isomer_b,
                "energy_a_hartree": e_a,
                "energy_b_hartree": e_b,
                "delta_energy_a_minus_b_hartree": e_a - e_b,
                "delta_energy_a_minus_b_kcal_mol": (e_a - e_b) * AU_TO_KCAL,
                "norm_full_a": row_a["norm_full"],
                "norm_full_b": row_b["norm_full"],
                "norm_projected_a": row_a["norm_projected"],
                "norm_projected_b": row_b["norm_projected"],
                "norm_removed_a": removed_a,
                "norm_removed_b": removed_b,
                "delta_norm_removed_a_minus_b": removed_a - removed_b,
                "lambda_x": row_a["lambda_x"],
                "lambda_y": row_a["lambda_y"],
                "lambda_z": row_a["lambda_z"],
            }
        )

    return pair_rows


def maybe_plot_pairwise(pair_rows: list[dict[str, object]], output_path: Path) -> None:
    if not pair_rows:
        return

    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("matplotlib is not available; skipping plot generation", file=sys.stderr)
        return

    phis = np.array([float(row["phi"]) for row in pair_rows])
    delta_energy = np.array([float(row["delta_energy_a_minus_b_kcal_mol"]) for row in pair_rows])
    removed_a = np.array([float(row["norm_removed_a"]) for row in pair_rows])
    removed_b = np.array([float(row["norm_removed_b"]) for row in pair_rows])
    isomer_a = str(pair_rows[0]["isomer_a"])
    isomer_b = str(pair_rows[0]["isomer_b"])
    theta = float(pair_rows[0]["theta"])

    fig, (ax_energy, ax_grad) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)

    ax_energy.plot(phis, delta_energy, marker="o", color="black")
    ax_energy.axhline(0.0, color="0.5", lw=1)
    ax_energy.set_ylabel(f"E({isomer_a}) - E({isomer_b})\n(kcal/mol)")
    ax_energy.set_title(f"theta = {theta:g}: energy difference and removed-gradient norm")
    ax_energy.grid(True, linestyle="--", alpha=0.35)

    ax_grad.plot(phis, removed_a, marker="o", label=f"{isomer_a}: |g_full - g_proj|")
    ax_grad.plot(phis, removed_b, marker="s", label=f"{isomer_b}: |g_full - g_proj|")
    ax_grad.set_xlabel("phi (deg.)")
    ax_grad.set_ylabel("removed-gradient norm")
    ax_grad.legend()
    ax_grad.grid(True, linestyle="--", alpha=0.35)

    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def configure_psi4(args) -> None:
    import psi4

    psi4.set_memory(args.memory)
    psi4.core.set_num_threads(args.num_threads)
    if args.psi4_output:
        psi4.core.set_output_file(str(args.psi4_output), False)
    else:
        psi4.core.be_quiet()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute full/projected CQED-DFT gradients for orientation rings."
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="runs_grid_no_freq folder")
    parser.add_argument("--thetas", default="72", help="Comma list or start:stop:step, e.g. 72 or 63:81:9")
    parser.add_argument("--phis", default="180:270:18", help="Comma list or start:stop:step")
    parser.add_argument("--isomers", nargs="+", default=list(DEFAULT_ISOMERS), help="Usually two isomers")
    parser.add_argument("--output-prefix", default="projected_gradient_ring", help="Output filename prefix")
    parser.add_argument("--dry-run", action="store_true", help="Validate folder parsing without Psi4 calculations")
    parser.add_argument("--plot", action="store_true", help="Write a quick pairwise diagnostic PNG")

    parser.add_argument("--memory", default="4 GB")
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument("--basis", default=DEFAULT_PSI4_OPTIONS["basis"])
    parser.add_argument("--scf-type", default=DEFAULT_PSI4_OPTIONS["scf_type"])
    parser.add_argument("--e-convergence", type=float, default=DEFAULT_PSI4_OPTIONS["e_convergence"])
    parser.add_argument("--d-convergence", type=float, default=DEFAULT_PSI4_OPTIONS["d_convergence"])
    parser.add_argument("--omega", type=float, default=0.06615)
    parser.add_argument("--reference", default="rks")
    parser.add_argument("--functional", default="wb97x")
    parser.add_argument("--charge", type=int, default=1)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--dispersion-policy", default="none")
    parser.add_argument("--density-fitting", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--quiet", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--psi4-output", type=Path, default=None, help="Optional Psi4 output file")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = args.root.resolve()
    thetas = parse_numeric_list(args.thetas)
    phis = parse_numeric_list(args.phis)
    cases = discover_cases(root, args.isomers, thetas, phis)
    warn_if_lambda_mismatch(cases)

    print(f"Discovered {len(cases)} cases in {root}")
    for case in cases:
        print(
            f"  {case.isomer:>5s} theta={case.theta:g} phi={case.phi:g} "
            f"lambda=({case.lambda_vector[0]: .8f}, {case.lambda_vector[1]: .8f}, {case.lambda_vector[2]: .8f})"
        )

    if args.dry_run:
        print("Dry run complete; no gradients were computed.")
        return

    configure_psi4(args)

    rows = []
    for index, case in enumerate(cases, start=1):
        print(f"[{index}/{len(cases)}] Computing {case.isomer} theta={case.theta:g} phi={case.phi:g}")
        rows.append(compute_case(case, args))

    raw_path = root / f"{args.output_prefix}_raw.csv"
    raw_fields = [
        "theta",
        "phi",
        "isomer",
        "energy_full_hartree",
        "energy_projected_hartree",
        "delta_energy_projected_minus_full_hartree",
        "norm_full",
        "norm_projected",
        "norm_removed",
        "lambda_x",
        "lambda_y",
        "lambda_z",
        "geometry_path",
    ]
    write_csv(rows, raw_path, raw_fields)
    print(f"Wrote {raw_path}")

    if len(args.isomers) >= 2:
        isomer_a, isomer_b = args.isomers[:2]
        pair_rows = build_pairwise_rows(rows, isomer_a, isomer_b)
        pair_path = root / f"{args.output_prefix}_{isomer_a}_{isomer_b}_pairwise.csv"
        pair_fields = [
            "theta",
            "phi",
            "isomer_a",
            "isomer_b",
            "energy_a_hartree",
            "energy_b_hartree",
            "delta_energy_a_minus_b_hartree",
            "delta_energy_a_minus_b_kcal_mol",
            "norm_full_a",
            "norm_full_b",
            "norm_projected_a",
            "norm_projected_b",
            "norm_removed_a",
            "norm_removed_b",
            "delta_norm_removed_a_minus_b",
            "lambda_x",
            "lambda_y",
            "lambda_z",
        ]
        write_csv(pair_rows, pair_path, pair_fields)
        print(f"Wrote {pair_path}")

        if args.plot:
            plot_path = root / f"{args.output_prefix}_{isomer_a}_{isomer_b}_diagnostic.png"
            maybe_plot_pairwise(pair_rows, plot_path)
            print(f"Wrote {plot_path}")


if __name__ == "__main__":
    main()
