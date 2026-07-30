#!/usr/bin/env python3
"""Validate completed QED-DFT optimization folders and optional spot recomputes.

Place this script in a campaign directory containing folders named like
``para_th18_ph108``. The fast checks do not import Psi4/CQED. Expensive
energy/gradient recomputes are only run for cells selected by command-line
options.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


CELL_RE = re.compile(r"^(para|meta|ortho)_th([-+]?\d+(?:\.\d+)?)_ph([-+]?\d+(?:\.\d+)?)$")
G_RE = re.compile(r"\|\s*g_proj\s*\|\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)")
E_RE = re.compile(r"\bE\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)")


@dataclass
class CellReport:
    cell_id: str
    path: Path
    status: str = "PASS"
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    recomputed: bool = False

    def fail(self, message: str) -> None:
        self.status = "FAIL"
        self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)


def load_json(path: Path) -> Any:
    with path.open() as handle:
        return json.load(handle)


def parse_xyz(path: Path) -> tuple[list[tuple[str, np.ndarray]], float | None, float | None]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError("XYZ file has fewer than two lines")

    try:
        n_atoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"invalid XYZ atom count: {lines[0]!r}") from exc

    if len(lines) < n_atoms + 2:
        raise ValueError(f"XYZ file has {len(lines) - 2} coordinate lines, expected {n_atoms}")

    comment = lines[1].strip()
    g_match = G_RE.search(comment)
    e_match = E_RE.search(comment)
    g_proj = float(g_match.group(1)) if g_match else None
    energy = float(e_match.group(1)) if e_match else None

    atoms: list[tuple[str, np.ndarray]] = []
    for idx, line in enumerate(lines[2 : 2 + n_atoms], start=1):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"coordinate line {idx} has fewer than four fields")
        atoms.append((parts[0], np.array([float(parts[1]), float(parts[2]), float(parts[3])])))

    return atoms, g_proj, energy


def field_vector(theta_deg: float, phi_deg: float) -> np.ndarray:
    theta = math.radians(theta_deg)
    phi = math.radians(phi_deg)
    return np.array(
        [
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(theta),
        ]
    )


def assert_close(
    report: CellReport,
    label: str,
    actual: Any,
    expected: Any,
    *,
    atol: float,
) -> None:
    if not np.allclose(np.array(actual, dtype=float), np.array(expected, dtype=float), atol=atol, rtol=0):
        report.fail(f"{label} mismatch: actual={actual}, expected={expected}, atol={atol}")


def validate_layout(cell_dir: Path, args: argparse.Namespace) -> tuple[CellReport, dict[str, Any] | None]:
    report = CellReport(cell_id=cell_dir.name, path=cell_dir)
    required_names = ["DONE", "cell.json", "opt_status.json", "opt_traj.xyz", "optimized.xyz"]
    missing = [name for name in required_names if not (cell_dir / name).is_file()]

    ccsd_files = sorted(cell_dir.glob("*_qed_ccsd_input.json"))
    if not ccsd_files:
        missing.append("*_qed_ccsd_input.json")
    elif len(ccsd_files) > 1:
        report.fail(f"found multiple *_qed_ccsd_input.json files: {[p.name for p in ccsd_files]}")

    if missing:
        message = f"missing required files: {', '.join(missing)}"
        if args.skip_incomplete:
            report.status = "SKIP"
            report.warnings.append(message)
            return report, None
        report.fail(message)
        return report, None

    try:
        cell = load_json(cell_dir / "cell.json")
        opt_status = load_json(cell_dir / "opt_status.json")
        ccsd = load_json(ccsd_files[0])
        xyz_atoms, xyz_g_proj, xyz_energy = parse_xyz(cell_dir / "optimized.xyz")
    except Exception as exc:
        report.fail(f"failed to read required data: {exc}")
        return report, None

    data = {
        "cell": cell,
        "opt_status": opt_status,
        "ccsd": ccsd,
        "xyz_atoms": xyz_atoms,
        "xyz_g_proj": xyz_g_proj,
        "xyz_energy": xyz_energy,
        "ccsd_path": ccsd_files[0],
    }
    return report, data


def validate_metadata(report: CellReport, data: dict[str, Any], args: argparse.Namespace) -> None:
    cell = data["cell"]
    opt_status = data["opt_status"]
    ccsd = data["ccsd"]
    xyz_atoms = data["xyz_atoms"]

    if cell.get("id") != report.cell_id:
        report.fail(f"cell.json id {cell.get('id')!r} does not match folder name")
    if opt_status.get("cell_id") != report.cell_id:
        report.fail(f"opt_status cell_id {opt_status.get('cell_id')!r} does not match folder name")

    match = CELL_RE.match(report.cell_id)
    if not match:
        report.fail("folder name does not match expected {isomer}_th{theta}_ph{phi} pattern")
    else:
        isomer, theta_text, phi_text = match.groups()
        if cell.get("isomer") != isomer:
            report.fail(f"cell.json isomer {cell.get('isomer')!r} does not match folder name")
        assert_close(report, "cell theta", cell.get("theta"), float(theta_text), atol=args.scalar_tol)
        assert_close(report, "cell phi", cell.get("phi"), float(phi_text), atol=args.scalar_tol)

    try:
        theta = float(cell.get("theta"))
        phi = float(cell.get("phi"))
        magnitude = float(cell.get("magnitude"))
    except (TypeError, ValueError) as exc:
        report.fail(f"cell.json theta/phi/magnitude is not numeric: {exc}")
        return
    expected_polvec = field_vector(theta, phi)
    expected_lambda = magnitude * expected_polvec

    assert_close(report, "cell lambda_vector", cell.get("lambda_vector"), expected_lambda, atol=args.vector_tol)
    assert_close(report, "opt_status lambda_vector", opt_status.get("lambda_vector"), expected_lambda, atol=args.vector_tol)

    if not opt_status.get("converged", False):
        report.fail("opt_status reports converged=false")
    final_gnorm = opt_status.get("final_gnorm")
    conv_threshold = opt_status.get("conv_threshold")
    if final_gnorm is not None and conv_threshold is not None and final_gnorm > conv_threshold:
        report.fail(f"final_gnorm {final_gnorm} exceeds conv_threshold {conv_threshold}")

    xyz_g_proj = data["xyz_g_proj"]
    xyz_energy = data["xyz_energy"]
    if xyz_g_proj is None:
        report.fail("optimized.xyz comment is missing |g_proj|")
    elif final_gnorm is not None:
        assert_close(report, "optimized.xyz |g_proj| vs opt_status final_gnorm", xyz_g_proj, final_gnorm, atol=args.gnorm_tol)
    if xyz_energy is None:
        report.fail("optimized.xyz comment is missing energy")
    elif opt_status.get("final_energy_hartree") is not None:
        assert_close(
            report,
            "optimized.xyz energy vs opt_status final_energy_hartree",
            xyz_energy,
            opt_status["final_energy_hartree"],
            atol=args.energy_tol,
        )

    geometry = ccsd.get("geometry", {})
    if geometry.get("units") != "angstrom":
        report.fail(f"CCSD geometry units are {geometry.get('units')!r}, expected 'angstrom'")
    if geometry.get("noorient") is not True:
        report.fail("CCSD geometry noorient is not true")

    ccsd_coords = geometry.get("coordinates", [])
    if len(ccsd_coords) != len(xyz_atoms):
        report.fail(f"CCSD coordinate count {len(ccsd_coords)} != optimized.xyz atom count {len(xyz_atoms)}")
    else:
        for idx, (ccsd_line, (xyz_elem, xyz_coord)) in enumerate(zip(ccsd_coords, xyz_atoms), start=1):
            parts = ccsd_line.split()
            if len(parts) < 4:
                report.fail(f"CCSD coordinate line {idx} has fewer than four fields")
                continue
            if parts[0] != xyz_elem:
                report.fail(f"atom {idx} element mismatch: CCSD={parts[0]}, XYZ={xyz_elem}")
                continue
            coord = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            if not np.allclose(coord, xyz_coord, atol=args.coord_tol, rtol=0):
                report.fail(f"atom {idx} coordinate mismatch beyond {args.coord_tol}")

    basis = ccsd.get("basis", {}).get("basisset")
    if basis != args.expected_basis:
        report.fail(f"CCSD basis {basis!r} != {args.expected_basis!r}")

    scf = ccsd.get("SCF", {})
    if scf.get("charge") != args.expected_charge:
        report.fail(f"CCSD charge {scf.get('charge')!r} != {args.expected_charge}")
    if scf.get("multiplicity") != args.expected_multiplicity:
        report.fail(f"CCSD multiplicity {scf.get('multiplicity')!r} != {args.expected_multiplicity}")
    qed_omegas = scf.get("qed_omegas", [])
    if len(qed_omegas) != 1:
        report.fail(f"CCSD qed_omegas has length {len(qed_omegas)}, expected 1")
    else:
        assert_close(report, "CCSD qed_omegas", qed_omegas[0], args.expected_omega, atol=args.scalar_tol)

    qed_lambdas = scf.get("qed_lambdas", [])
    if len(qed_lambdas) != 1:
        report.fail(f"CCSD qed_lambdas has length {len(qed_lambdas)}, expected 1")
    else:
        assert_close(report, "CCSD qed_lambdas", qed_lambdas[0], magnitude, atol=args.scalar_tol)

    polvecs = scf.get("qed_polvecs", [])
    if len(polvecs) != 1:
        report.fail(f"CCSD qed_polvecs has length {len(polvecs)}, expected 1")
    else:
        assert_close(report, "CCSD qed_polvecs", polvecs[0], expected_polvec, atol=args.vector_tol)

    task = ccsd.get("TASK", {})
    if task.get("ccsd") is not True:
        report.fail("CCSD TASK.ccsd is not true")
    if task.get("scf") is not False:
        report.fail("CCSD TASK.scf is not false")


def build_geometry_string(xyz_atoms: list[tuple[str, np.ndarray]], charge: int, multiplicity: int) -> str:
    lines = [f"{elem:4s} {coord[0]:20.12f} {coord[1]:20.12f} {coord[2]:20.12f}" for elem, coord in xyz_atoms]
    return f"""{charge} {multiplicity}
{chr(10).join(lines)}
units angstrom
no_reorient
no_com
symmetry c1
"""


def recompute_cell(report: CellReport, data: dict[str, Any], args: argparse.Namespace) -> None:
    try:
        import psi4
        from cqed_scf.calculator import CQEDCalculator
    except Exception as exc:
        report.fail(f"could not import Psi4/CQED recompute dependencies: {exc}")
        return

    cell = data["cell"]
    xyz_energy = data["xyz_energy"]
    xyz_g_proj = data["xyz_g_proj"]
    if xyz_energy is None or xyz_g_proj is None:
        report.fail("cannot recompute without energy and |g_proj| in optimized.xyz comment")
        return

    psi4_options = {
        "basis": args.expected_basis,
        "reference": "rks",
        "scf_type": "df",
        "e_convergence": 1e-9,
        "d_convergence": 1e-9,
        "dft_radial_points": 90,
        "dft_spherical_points": 590,
        "dft_pruning_scheme": "none",
    }

    lambda_vector = np.array(cell["lambda_vector"], dtype=float)
    geometry_string = build_geometry_string(data["xyz_atoms"], args.expected_charge, args.expected_multiplicity)
    output_path = report.path / args.recompute_output

    old_cwd = Path.cwd()
    try:
        os.chdir(report.path)
        psi4.core.set_output_file(str(output_path), False)
        calc = CQEDCalculator(
            lambda_vector=lambda_vector,
            psi4_options=psi4_options,
            omega=args.expected_omega,
            density_fitting=True,
            charge=args.expected_charge,
            multiplicity=args.expected_multiplicity,
            functional=args.functional,
            reference="rks",
            dispersion_policy="none",
            debug=False,
        )
        energy, gradient, _ = calc.energy_and_projected_gradient(geometry_string)
    except Exception as exc:
        report.fail(f"recompute failed: {exc}")
        return
    finally:
        os.chdir(old_cwd)

    g_proj = float(np.linalg.norm(gradient))
    report.recomputed = True
    assert_close(report, "recomputed energy", energy, xyz_energy, atol=args.recompute_energy_tol)
    assert_close(report, "recomputed |g_proj|", g_proj, xyz_g_proj, atol=args.recompute_gnorm_tol)


def discover_cells(root: Path) -> list[Path]:
    return sorted(path for path in root.iterdir() if path.is_dir() and CELL_RE.match(path.name))


def select_recomputes(passed: list[CellReport], args: argparse.Namespace) -> set[str]:
    selected: set[str] = set()
    by_id = {report.cell_id: report for report in passed}

    if args.recompute_all:
        selected.update(by_id)

    for item in args.recompute_cell:
        for cell_id in item.split(","):
            cell_id = cell_id.strip()
            if cell_id:
                selected.add(cell_id)

    rng = random.Random(args.seed)
    remaining = [report.cell_id for report in passed if report.cell_id not in selected]
    if args.recompute_probability > 0:
        selected.update(cell_id for cell_id in remaining if rng.random() < args.recompute_probability)

    remaining = [report.cell_id for report in passed if report.cell_id not in selected]
    if args.recompute_count > 0:
        selected.update(rng.sample(remaining, min(args.recompute_count, len(remaining))))

    return selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", default=".", type=Path, help="campaign directory to validate")
    parser.add_argument("--skip-incomplete", action="store_true", help="skip folders missing required completion files")
    parser.add_argument("--json-report", type=Path, help="optional path for a machine-readable report")

    parser.add_argument("--expected-basis", default="6-311G*")
    parser.add_argument("--expected-charge", type=int, default=1)
    parser.add_argument("--expected-multiplicity", type=int, default=1)
    parser.add_argument("--expected-omega", type=float, default=0.06615)
    parser.add_argument("--functional", default="wb97x")

    parser.add_argument("--coord-tol", type=float, default=1e-10)
    parser.add_argument("--vector-tol", type=float, default=1e-12)
    parser.add_argument("--scalar-tol", type=float, default=1e-12)
    parser.add_argument("--energy-tol", type=float, default=5e-10)
    parser.add_argument("--gnorm-tol", type=float, default=5e-7)

    parser.add_argument("--recompute-cell", action="append", default=[], help="cell id to recompute; may be repeated or comma-separated")
    parser.add_argument("--recompute-all", action="store_true", help="recompute every cell that passes fast checks")
    parser.add_argument("--recompute-count", type=int, default=0, help="randomly recompute this many fast-passing cells")
    parser.add_argument("--recompute-probability", type=float, default=0.0, help="randomly recompute each fast-passing cell with this probability")
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--recompute-energy-tol", type=float, default=1e-8)
    parser.add_argument("--recompute-gnorm-tol", type=float, default=1e-6)
    parser.add_argument("--recompute-output", default="validation_recompute.out")

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        print(f"ERROR: campaign root does not exist: {root}", file=sys.stderr)
        return 2

    cell_dirs = discover_cells(root)
    if not cell_dirs:
        print(f"ERROR: found no orientation folders below {root}", file=sys.stderr)
        return 2

    reports: list[CellReport] = []
    data_by_id: dict[str, dict[str, Any]] = {}
    for cell_dir in cell_dirs:
        report, data = validate_layout(cell_dir, args)
        if data is not None:
            validate_metadata(report, data, args)
            data_by_id[report.cell_id] = data
        reports.append(report)

    fast_passed = [report for report in reports if report.status == "PASS"]
    selected = select_recomputes(fast_passed, args)
    unknown = sorted(cell_id for cell_id in selected if cell_id not in data_by_id)
    for cell_id in unknown:
        print(f"WARNING: requested recompute cell was not found or did not pass fast checks: {cell_id}")
        selected.remove(cell_id)

    if selected:
        print(f"Running expensive recomputes for {len(selected)} cell(s): {', '.join(sorted(selected))}")
    for report in reports:
        if report.cell_id in selected and report.status == "PASS":
            recompute_cell(report, data_by_id[report.cell_id], args)

    counts = {"PASS": 0, "FAIL": 0, "SKIP": 0}
    for report in reports:
        counts[report.status] = counts.get(report.status, 0) + 1

    print(f"Validated {len(reports)} orientation folder(s) below {root}")
    print(f"PASS: {counts.get('PASS', 0)}  FAIL: {counts.get('FAIL', 0)}  SKIP: {counts.get('SKIP', 0)}")
    recomputed = [report.cell_id for report in reports if report.recomputed]
    if recomputed:
        print(f"Recomputed: {', '.join(sorted(recomputed))}")

    for report in reports:
        if report.status == "FAIL":
            print(f"\nFAIL {report.cell_id}")
            for error in report.errors:
                print(f"  - {error}")
        elif report.warnings and not args.skip_incomplete:
            print(f"\nWARN {report.cell_id}")
            for warning in report.warnings:
                print(f"  - {warning}")

    if args.json_report:
        payload = [
            {
                "cell_id": report.cell_id,
                "path": str(report.path),
                "status": report.status,
                "errors": report.errors,
                "warnings": report.warnings,
                "recomputed": report.recomputed,
            }
            for report in reports
        ]
        args.json_report.write_text(json.dumps(payload, indent=2) + "\n")

    return 1 if counts.get("FAIL", 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())
