"""
ccsd_input.py

Builds a QED-CCSD input JSON (the format used to launch the downstream
frequency-independent QED-CCSD/EOM-CCSD calculation) from a converged
constrained-geometry-optimization cell.

This is a separate module from cavity_common.py on purpose: cavity_common.py
is the single source of truth for the QED-DFT *optimization* physics (what
grid_campaign_no_freq.py / campaign.py / trajectory_campaign.py actually run);
this module only formats the *output* of that optimization into an input file
for a different downstream program (QED-CCSD), at a different, fixed lambda
magnitude convention. Keeping it separate means a change to the DFT opt
options never accidentally perturbs the CCSD input template, and vice versa.

Only two things vary per cell:
  - geometry.coordinates      (the converged geometry for this cell)
  - SCF.qed_polvecs            (the unit field-direction vector for this
                                cell's theta, phi)

Everything else is fixed, taken either literally from the example input file
Jay provided, or reused from cavity_common.py where the same physical
constant already has a single source of truth (omega, basis, charge,
multiplicity).
"""

import json

# cavity_common is imported lazily (inside the functions that need it), not
# at module level, matching campaign.py's own convention -- cavity_common
# imports psi4/ase/cqed_scf eagerly, and callers that only want to build a
# manifest or dry-run a campaign (no actual optimization) shouldn't need
# those installed just to import this module.

# The fixed cavity coupling magnitude used throughout this campaign. All
# QED-CCSD input files generated from these optimizations report this same
# qed_lambdas value, regardless of theta/phi (only the *direction*,
# qed_polvecs, varies from cell to cell).
QED_CCSD_LAMBDA_MAG = 0.1


def format_geometry_coordinates(symbols, coords_angstrom):
    """
    Render (symbols, coords_angstrom) as the list-of-strings format the
    QED-CCSD input expects, e.g. "C   -0.493854901900   1.213549282200   0.828747960700".
    """
    lines = []
    for sym, (x, y, z) in zip(symbols, coords_angstrom):
        lines.append(f"{sym:<2s}  {x:16.12f}  {y:16.12f}  {z:16.12f}")
    return lines


def qed_polvec_for(theta_deg, phi_deg):
    """
    Unit field-direction vector for (theta, phi), reusing cavity_common's
    single lambda-vector construction (magnitude=1.0 gives the bare unit
    direction -- the same helper that builds the actual optimization's
    lambda_vector, just unscaled).
    """
    import cavity_common as cc
    return cc.lambda_vector_for(theta_deg, phi_deg, magnitude=1.0)


def build_qed_ccsd_input(theta_deg, phi_deg, symbols, coords_angstrom):
    """
    Build the full QED-CCSD input dict for one converged cell.

    Parameters
    ----------
    theta_deg, phi_deg : float
        The cell's cavity-orientation angles (degrees).
    symbols : list[str]
        Atomic symbols, in the same order as coords_angstrom.
    coords_angstrom : (natom, 3) array-like
        The converged (optimized) Cartesian geometry, in Angstrom.

    Returns
    -------
    dict
        Ready to json.dump directly to a QED-CCSD input file.
    """
    import cavity_common as cc

    polvec = qed_polvec_for(theta_deg, phi_deg)

    return {
        "geometry": {
            "coordinates": format_geometry_coordinates(symbols, coords_angstrom),
            "units": "angstrom",
            "noorient": True,
        },
        "common": {
            "maxiter": 100,
            "debug": False,
        },
        "basis": {
            "basisset": cc.PSI4_OPTIONS["basis"],
        },
        "SCF": {
            "charge": cc.CHARGE,
            "multiplicity": cc.MULTIPLICITY,
            "lshift": 0,
            "tol_lindep": 1e-10,
            "conve": 1e-07,
            "convd": 1e-05,
            "diis_hist": 10,
            "writem": 10,
            "restart": False,
            "noscf": False,
            "scf_type": "restricted",
            "debug": False,
            "qed_omegas": [cc.OMEGA],
            "qed_volumes": [],
            "qed_lambdas": [QED_CCSD_LAMBDA_MAG],
            "qed_polvecs": [polvec],
            "direct_df": False,
        },
        "CD": {
            "debug": False,
            "diagtol": 1e-06,
        },
        "CC": {
            "threshold": 1e-06,
            "ndiis": 5,
            "ccsd_maxiter": 200,
            "tilesize": 40,
            "freeze": {
                "atomic": True,
            },
            "readt": False,
            "writet": True,
            "writet_iter": 5,
            "debug": False,
            "profile_ccsd": False,
            "balance_tiles": True,
        },
        "TASK": {
            "scf": False,
            "ccsd": True,
        },
    }


def write_qed_ccsd_input(path, theta_deg, phi_deg, symbols, coords_angstrom):
    """Build and write the QED-CCSD input JSON for one converged cell to `path`."""
    data = build_qed_ccsd_input(theta_deg, phi_deg, symbols, coords_angstrom)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    return data
