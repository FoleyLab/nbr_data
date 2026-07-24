#!/usr/bin/env python3
"""
Pauli-Fierz Hamiltonian lambda-magnitude scan for a fixed (theta, phi)
cavity-field direction, for the ortho, meta, and para nitro-isomers.

Companion script to pf_isomer_scan.py.  Where pf_isomer_scan.py scans over a
grid of (theta, phi) directions at a single fixed lambda magnitude,
this script fixes the (theta, phi) direction and instead scans over a list
of lambda magnitudes -- mirroring the scan performed by the CQED-SCF input
scripts in nitrobenzene_opt_and_freq/input_scripts (e.g.
input_ortho_meta_dir70_31.py and input_para_meta_dir65_78.py), but using the
Pauli-Fierz Hamiltonian parameterized by EOM-CCSD data instead of CQED-SCF.

All EOM-CCSD data parsing and PF Hamiltonian construction is reused directly
from pf_isomer_scan.py (imported as a module, not modified/duplicated).

Outputs (CSV files, written alongside this script):
  pqed_<Nel>_<Nph>_dir<theta>_<phi>_scan.csv          total energies, non-CS
  pqed_<Nel>_<Nph>_dir<theta>_<phi>_scan_CS.csv        total energies, CS
  pqed_<Nel>_<Nph>_dir<theta>_<phi>_scan_decomp.csv    energy decomposition, non-CS  (only if WRITE_DECOMPOSITION)
  pqed_<Nel>_<Nph>_dir<theta>_<phi>_scan_decomp_CS.csv energy decomposition, CS       (only if WRITE_DECOMPOSITION)

Only the configuration block below should normally need editing.
"""

import os
import sys
import csv
import numpy as np

# Make sure pf_isomer_scan.py (expected to live next to this file) is
# importable regardless of the current working directory.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import pf_isomer_scan as pfscan  # noqa: E402  (reuse data parsing / PF Hamiltonian builder)

# ============================================================================
#                              CONFIGURATION
# ============================================================================

# --- PF Hamiltonian size knobs ----------------------------------------------
NUM_ELECTRONIC_STATES = 49   # number of electronic states used to build H_PF
NUM_FOCK_STATES       = 3   # number of photonic Fock states

# --- Photon frequency (Hartree), fit from the EOM-CCSD data -----------------
OMEGA = 0.066148

# --- Fixed cavity-field direction (degrees) ---------------------------------
THETA = 70.0   # degrees from the z-axis
PHI   = 31.0   # degrees from the x-axis in the xy-plane

# --- Lambda magnitudes to scan over, same (theta, phi) direction for all ----
LAMBDA_MAGNITUDES = [0.02, 0.04, 0.06, 0.08, 0.10]

# --- Set True to also write the (E_el, E_ph, E_blc, E_dse) decomposition ---
# --- CSVs (one for non-CS, one for CS). Set False to skip them. ------------
WRITE_DECOMPOSITION = True

HARTREE_TO_KCAL_MOL = 627.5094740631

# --- Canonical isomer output order (independent of pfscan.ISOMERS order) ---
LABELS = ["Ortho", "Meta", "Para"]

# --- Output file names, derived from the config knobs above ----------------
_TAG = f"{NUM_ELECTRONIC_STATES}_{NUM_FOCK_STATES}_dir{int(THETA)}_{int(PHI)}"
TOTAL_ENERGY_CSV    = os.path.join(SCRIPT_DIR, f"pqed_{_TAG}_scan.csv")
TOTAL_ENERGY_CSV_CS = os.path.join(SCRIPT_DIR, f"pqed_{_TAG}_scan_CS.csv")
DECOMP_CSV          = os.path.join(SCRIPT_DIR, f"pqed_{_TAG}_scan_decomp.csv")
DECOMP_CSV_CS       = os.path.join(SCRIPT_DIR, f"pqed_{_TAG}_scan_decomp_CS.csv")


# ============================================================================
#                              DATA LOADING
# ============================================================================

def load_isomer_data():
    """Parse EOM-CCSD data for each isomer via pf_isomer_scan.py's parser."""
    isomer_data = {}
    for label, fname in pfscan.ISOMERS.items():
        path = os.path.join(SCRIPT_DIR, pfscan.DATA_DIR, fname)
        ref_E, corr_E, excit_E, electronic_dipole, total_dipoles = pfscan.parse_cq_h5_data(
            path, verbose=True
        )

        # symmetrize the transition-dipole matrix (bra/ket swap), as in pf_isomer_scan.py
        electronic_dipoles_sym = 0.5 * (
            electronic_dipole + np.transpose(electronic_dipole, axes=(1, 0, 2))
        )
        total_dipoles_sym = 0.5 * (
            total_dipoles + np.transpose(total_dipoles, axes=(1, 0, 2))
        )

        E_el = pfscan.build_electronic_energies(ref_E, corr_E, excit_E)

        if NUM_ELECTRONIC_STATES > len(E_el):
            raise ValueError(
                f"{label}: NUM_ELECTRONIC_STATES={NUM_ELECTRONIC_STATES} "
                f"exceeds available states ({len(E_el)})."
            )

        isomer_data[label] = {
            "E_el": E_el,
            "electronic_dipoles": electronic_dipoles_sym,
            "total_dipoles": total_dipoles_sym,
        }
    return isomer_data


# ============================================================================
#                              OUTPUT HELPERS
# ============================================================================

def _write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# ============================================================================
#                                 MAIN
# ============================================================================

def main():
    isomer_data = load_isomer_data()

    lambda_dir = pfscan.generate_lambda_vec_from_theta_and_phi(THETA, PHI)

    total_rows, total_cs_rows = [], []
    decomp_rows, decomp_cs_rows = [], []

    print(f"Running PF lambda-magnitude scan at theta={THETA}, phi={PHI}")
    print(
        f"  electronic states = {NUM_ELECTRONIC_STATES}, "
        f"Fock states = {NUM_FOCK_STATES}, omega = {OMEGA}"
    )
    print(f"  lambda magnitudes = {LAMBDA_MAGNITUDES}")

    for lam_mag in LAMBDA_MAGNITUDES:
        lam = lambda_dir * lam_mag
        Ex, Ey, Ez = lam

        totals, totals_cs = {}, {}
        decomp, decomp_cs = {}, {}

        for label in LABELS:
            d = isomer_data[label]
            (
                H_total, H_total_cs, H_el_full, H_ph_full,
                H_blc, H_dse, H_blc_cs, H_dse_cs,
            ) = pfscan.build_PF_Hamiltonian(
                dim_ph=NUM_FOCK_STATES,
                dim_el=NUM_ELECTRONIC_STATES,
                omega=OMEGA,
                lambda_vec=lam,
                e_el=d["E_el"],
                electronic_dipoles=d["electronic_dipoles"],
                total_dipoles=d["total_dipoles"],
            )

            # non-CS: diagonalize H_total (total-dipole H_blc/H_dse)
            eigs, C = np.linalg.eigh(H_total)
            psi0 = C[:, 0]

            # CS: diagonalize H_total_cs (electronic-dipole-only, shifted H_blc_cs/H_dse_cs)
            eigs_cs, C_cs = np.linalg.eigh(H_total_cs)
            psi0_cs = C_cs[:, 0]

            totals[label] = eigs[0]
            totals_cs[label] = eigs_cs[0]

            if WRITE_DECOMPOSITION:
                decomp[label] = (
                    np.real(psi0.conj().T @ H_el_full @ psi0),
                    np.real(psi0.conj().T @ H_ph_full @ psi0),
                    np.real(psi0.conj().T @ H_blc @ psi0),
                    np.real(psi0.conj().T @ H_dse @ psi0),
                )
                decomp_cs[label] = (
                    np.real(psi0_cs.conj().T @ H_el_full @ psi0_cs),
                    np.real(psi0_cs.conj().T @ H_ph_full @ psi0_cs),
                    np.real(psi0_cs.conj().T @ H_blc_cs @ psi0_cs),
                    np.real(psi0_cs.conj().T @ H_dse_cs @ psi0_cs),
                )

        # pairwise energy differences (kcal/mol)
        dE_ortho_meta = (totals["Ortho"] - totals["Meta"]) * HARTREE_TO_KCAL_MOL
        dE_para_meta = (totals["Para"] - totals["Meta"]) * HARTREE_TO_KCAL_MOL
        dE_ortho_para = (totals["Ortho"] - totals["Para"]) * HARTREE_TO_KCAL_MOL

        dE_ortho_meta_cs = (totals_cs["Ortho"] - totals_cs["Meta"]) * HARTREE_TO_KCAL_MOL
        dE_para_meta_cs = (totals_cs["Para"] - totals_cs["Meta"]) * HARTREE_TO_KCAL_MOL
        dE_ortho_para_cs = (totals_cs["Ortho"] - totals_cs["Para"]) * HARTREE_TO_KCAL_MOL

        row = {
            "theta": THETA,
            "phi": PHI,
            "Ex": Ex,
            "Ey": Ey,
            "Ez": Ez,
            "lambda_magnitude": lam_mag,
            "E_ortho_Hartrees": totals["Ortho"],
            "E_meta_Hartrees": totals["Meta"],
            "E_para_Hartrees": totals["Para"],
            "dE_ortho_meta_kcal/mol": dE_ortho_meta,
            "dE_para_meta_kcal/mol": dE_para_meta,
            "dE_ortho_para_kcal/mol": dE_ortho_para,
        }
        row_cs = {
            **row,
            "E_ortho_Hartrees": totals_cs["Ortho"],
            "E_meta_Hartrees": totals_cs["Meta"],
            "E_para_Hartrees": totals_cs["Para"],
            "dE_ortho_meta_kcal/mol": dE_ortho_meta_cs,
            "dE_para_meta_kcal/mol": dE_para_meta_cs,
            "dE_ortho_para_kcal/mol": dE_ortho_para_cs,
        }

        total_rows.append(row)
        total_cs_rows.append(row_cs)

        if WRITE_DECOMPOSITION:
            drow = {
                "theta": THETA, "phi": PHI, "Ex": Ex, "Ey": Ey, "Ez": Ez,
                "lambda_magnitude": lam_mag,
            }
            drow_cs = dict(drow)
            for label in LABELS:
                l = label.lower()
                e_el, e_ph, e_blc, e_dse = decomp[label]
                drow[f"{l}_E_el"] = e_el
                drow[f"{l}_E_ph"] = e_ph
                drow[f"{l}_E_blc"] = e_blc
                drow[f"{l}_E_dse"] = e_dse

                e_el_cs, e_ph_cs, e_blc_cs, e_dse_cs = decomp_cs[label]
                drow_cs[f"{l}_E_el"] = e_el_cs
                drow_cs[f"{l}_E_ph"] = e_ph_cs
                drow_cs[f"{l}_E_blc"] = e_blc_cs
                drow_cs[f"{l}_E_dse"] = e_dse_cs

            decomp_rows.append(drow)
            decomp_cs_rows.append(drow_cs)

        print(
            f"  lambda_magnitude={lam_mag}: "
            f"E_ortho={totals['Ortho']:.8f}  E_meta={totals['Meta']:.8f}  "
            f"E_para={totals['Para']:.8f}  (non-CS)"
        )

    # --- write total-energy CSVs -------------------------------------------
    total_fields = [
        "theta", "phi", "Ex", "Ey", "Ez", "lambda_magnitude",
        "E_ortho_Hartrees", "E_meta_Hartrees", "E_para_Hartrees",
        "dE_ortho_meta_kcal/mol", "dE_para_meta_kcal/mol", "dE_ortho_para_kcal/mol",
    ]
    _write_csv(TOTAL_ENERGY_CSV, total_fields, total_rows)
    _write_csv(TOTAL_ENERGY_CSV_CS, total_fields, total_cs_rows)

    print(f"\nWrote {TOTAL_ENERGY_CSV}")
    print(f"Wrote {TOTAL_ENERGY_CSV_CS}")

    # --- write decomposition CSVs (optional) -------------------------------
    if WRITE_DECOMPOSITION:
        decomp_fields = ["theta", "phi", "Ex", "Ey", "Ez", "lambda_magnitude"]
        for label in LABELS:
            l = label.lower()
            decomp_fields += [f"{l}_E_el", f"{l}_E_ph", f"{l}_E_blc", f"{l}_E_dse"]

        _write_csv(DECOMP_CSV, decomp_fields, decomp_rows)
        _write_csv(DECOMP_CSV_CS, decomp_fields, decomp_cs_rows)

        print(f"Wrote {DECOMP_CSV}")
        print(f"Wrote {DECOMP_CSV_CS}")


if __name__ == "__main__":
    main()
