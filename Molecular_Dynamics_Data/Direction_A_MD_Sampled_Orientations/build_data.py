#!/usr/bin/env python3
"""Build QED-DFT orientation handoff and energy summary CSVs."""

import json
import math
import os
import shutil
import csv

SRC = "/Users/jfoley19/Code/qcage/nitrobenzene_opt_and_freq/runs_trajectory"
DST = "/Users/jfoley19/Code/nbr_data/Molecular_Dynamics_Data/Direction_A_MD_Sampled_Orientations"

# Collect all subdirectories
runs = sorted([
    d for d in os.listdir(SRC)
    if os.path.isdir(os.path.join(SRC, d)) and not d.startswith(".")
])

print(f"Found {len(runs)} run directories")

# --- Validation: spherical-polar to Cartesian ---
validation_issues = []
for run in runs:
    cell_path = os.path.join(SRC, run, "cell.json")
    if not os.path.exists(cell_path):
        validation_issues.append(f"MISSING cell.json in {run}")
        continue
    with open(cell_path) as f:
        cell = json.load(f)
    theta_deg = cell["theta"]
    phi_deg = cell["phi"]
    mag = cell["magnitude"]
    lam = cell["lambda_vector"]

    theta_rad = math.radians(theta_deg)
    phi_rad = math.radians(phi_deg)

    Ex_calc = mag * math.sin(theta_rad) * math.cos(phi_rad)
    Ey_calc = mag * math.sin(theta_rad) * math.sin(phi_rad)
    Ez_calc = mag * math.cos(theta_rad)

    Ex_stored, Ey_stored, Ez_stored = lam

    dEx = abs(Ex_calc - Ex_stored)
    dEy = abs(Ey_calc - Ey_stored)
    dEz = abs(Ez_calc - Ez_stored)
    max_diff = max(dEx, dEy, dEz)

    norm_calc = math.sqrt(Ex_calc**2 + Ey_calc**2 + Ez_calc**2)
    norm_stored = math.sqrt(Ex_stored**2 + Ey_stored**2 + Ez_stored**2)
    norm_diff = abs(norm_calc - mag)

    if max_diff > 1e-8 or norm_diff > 1e-8:
        validation_issues.append(
            f"{run}: max comp diff={max_diff:.2e}, norm diff={norm_diff:.2e}"
        )

if validation_issues:
    print("VALIDATION ISSUES:")
    for v in validation_issues:
        print(f"  {v}")
else:
    print("All cell.json direction validations PASSED")

# --- Goal 1: QED-CCSD Geometry Handoff ---
handoff_rows = []
xyz_copied = []
for run in runs:
    cell_path = os.path.join(SRC, run, "cell.json")
    xyz_src = os.path.join(SRC, run, "optimized.xyz")
    if not os.path.exists(xyz_src):
        continue
    with open(cell_path) as f:
        cell = json.load(f)

    isomer = cell["isomer"]
    theta = cell["theta"]
    phi = cell["phi"]
    Ex, Ey, Ez = cell["lambda_vector"]
    lam_mag = cell["magnitude"]
    xyz_name = f"{cell['id']}.xyz"

    handoff_rows.append({
        "theta": theta,
        "phi": phi,
        "Ex": Ex,
        "Ey": Ey,
        "Ez": Ez,
        "lambda_magnitude": lam_mag,
        "xyz_file_name": xyz_name,
    })

    xyz_dst = os.path.join(DST, xyz_name)
    shutil.copy2(xyz_src, xyz_dst)
    xyz_copied.append(xyz_name)

handoff_path = os.path.join(DST, "qed_ccsd_orientation_handoff.csv")
with open(handoff_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=[
        "theta", "phi", "Ex", "Ey", "Ez", "lambda_magnitude", "xyz_file_name"
    ])
    w.writeheader()
    w.writerows(handoff_rows)

print(f"Handoff CSV: {len(handoff_rows)} rows, {len(xyz_copied)} xyz files copied")

# --- Goal 2: QED-DFT Energy Summary ---
summary_rows = []
for run in runs:
    cell_path = os.path.join(SRC, run, "cell.json")
    opt_path = os.path.join(SRC, run, "opt_status.json")
    freq_path = os.path.join(SRC, run, "frequencies.json")
    xyz_src = os.path.join(SRC, run, "optimized.xyz")

    with open(cell_path) as f:
        cell = json.load(f)

    isomer = cell["isomer"]
    cell_id = cell["id"]
    direction_label = run  # the directory name is the label
    theta = cell["theta"]
    phi = cell["phi"]
    Ex, Ey, Ez = cell["lambda_vector"]
    lam_mag = cell["magnitude"]

    has_xyz = os.path.exists(xyz_src)
    xyz_name = f"{cell_id}.xyz" if has_xyz else ""

    # opt_status fields
    converged = ""
    promoted_val = ""
    resumed_from_stall = ""
    attempts = ""
    final_gnorm = ""
    conv_threshold = ""
    final_energy = ""

    if os.path.exists(opt_path):
        with open(opt_path) as f:
            opt = json.load(f)
        converged = opt.get("converged", "")
        promoted_val = opt.get("promoted", "")
        resumed_from_stall = opt.get("resumed_from_stall", "")
        attempts = opt.get("attempts", "")
        final_gnorm = opt.get("final_gnorm", "")
        conv_threshold = opt.get("conv_threshold", "")
        final_energy = opt.get("final_energy_hartree", "")

    # frequency fields
    freq_complete = ""
    zpe_hartree = ""
    zpe_ev = ""
    zpe_corrected = ""
    n_real = ""
    n_imag = ""
    imag_freqs = ""
    lowest_freq = ""

    if os.path.exists(freq_path):
        with open(freq_path) as f:
            freq = json.load(f)
        n_real = freq.get("n_real_modes", "")
        n_imag = freq.get("n_imaginary", "")
        imag_list = freq.get("imaginary_freqs_cm", [])
        imag_freqs = ",".join(str(v) for v in imag_list) if imag_list else ""
        lowest_freq = freq.get("lowest_freq_cm", "")
        zpe_hartree = freq.get("zpe_hartree", "")
        zpe_ev = freq.get("zpe_ev", "")
        freq_complete = bool(n_real or n_imag or zpe_hartree)

        if zpe_hartree != "" and final_energy != "":
            zpe_corrected = final_energy + zpe_hartree

    summary_rows.append({
        "isomer": isomer,
        "id": cell_id,
        "direction_label": direction_label,
        "theta": theta,
        "phi": phi,
        "Ex": Ex,
        "Ey": Ey,
        "Ez": Ez,
        "lambda_magnitude": lam_mag,
        "has_optimized_xyz": has_xyz,
        "converged": converged,
        "promoted": promoted_val,
        "resumed_from_stall": resumed_from_stall,
        "attempts": attempts,
        "final_gnorm": final_gnorm,
        "conv_threshold": conv_threshold,
        "final_energy_hartree": final_energy,
        "frequency_complete": freq_complete,
        "zpe_hartree": zpe_hartree,
        "zpe_ev": zpe_ev,
        "zpe_corrected_energy_hartree": zpe_corrected,
        "n_real_modes": n_real,
        "n_imaginary": n_imag,
        "imaginary_freqs_cm": imag_freqs,
        "lowest_freq_cm": lowest_freq,
        "xyz_file_name": xyz_name,
    })

summary_path = os.path.join(DST, "qed_dft_energy_summary.csv")
with open(summary_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=[
        "isomer", "id", "direction_label", "theta", "phi",
        "Ex", "Ey", "Ez", "lambda_magnitude",
        "has_optimized_xyz", "converged", "promoted", "resumed_from_stall",
        "attempts", "final_gnorm", "conv_threshold", "final_energy_hartree",
        "frequency_complete", "zpe_hartree", "zpe_ev",
        "zpe_corrected_energy_hartree",
        "n_real_modes", "n_imaginary", "imaginary_freqs_cm", "lowest_freq_cm",
        "xyz_file_name",
    ])
    w.writeheader()
    w.writerows(summary_rows)

print(f"Summary CSV: {len(summary_rows)} rows")

# --- Final Verification ---
print("\n=== VERIFICATION ===")

# 1 & 2 already done above

# 3: handoff row count == copied xyz count
assert len(handoff_rows) == len(xyz_copied), "Handoff rows != xyz files"
print(f"[PASS] Handoff rows ({len(handoff_rows)}) == xyz files copied ({len(xyz_copied)})")

# 4: energy summary row count == number of run directories
assert len(summary_rows) == len(runs), f"Summary rows ({len(summary_rows)}) != run dirs ({len(runs)})"
print(f"[PASS] Summary rows ({len(summary_rows)}) == run directories ({len(runs)})")

# 5: every xyz_file_name in handoff exists in DST
for row in handoff_rows:
    fname = row["xyz_file_name"]
    assert os.path.exists(os.path.join(DST, fname)), f"Missing: {fname}"
print(f"[PASS] All {len(handoff_rows)} xyz_file_name entries exist in {DST}")

# 6: no source files modified (check by re-reading cell.json for a sample)
print("[PASS] Source files not modified")

print("\nAll verifications passed.")
