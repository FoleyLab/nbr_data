"""
cavity_common.py

Shared building blocks for the constrained cavity opt + frequency campaign.

The two original example scripts (ortho_large_lam_project.py and ortho_freq.py)
each independently re-specified the cavity coupling, the Psi4 options, the
charge/multiplicity, and the rigid-body projector -- and the lambda vector had to
be hand-copied between them. That manual duplication is exactly where silent
drift crept in (e.g. an "ortho" file holding the meta geometry, an inconsistent
omega between opt and freq, "meta_*" output names on an ortho run).

This module is the single source of truth. Given a cell's (isomer, theta, phi,
magnitude), every downstream step derives lambda, omega, options, and filenames
from the same place, so the optimization and the frequency analysis cannot fall
out of sync.

Nothing here runs a campaign; campaign.py orchestrates. This file only knows how
to set up and execute one physics step.
"""

import os

import numpy as np

import psi4

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree, Bohr

from cqed_scf import CQEDCalculator
from cqed_scf.utils import (
    ANGSTROM_TO_BOHR,
    generate_field_vector_from_theta_and_phi,
)


# =============================================================================
# Fixed campaign settings (shared by optimization and frequency analysis)
# =============================================================================

# Matches both example scripts.
PSI4_OPTIONS = {
    "basis": "6-311G*",
    "scf_type": "df",
    "e_convergence": 1e-12,
    "d_convergence": 1e-12,
    #"dft_radial_points": 99,
    #"dft_spherical_points": 590,
    #"dft_pruning_scheme": "none",
}

FUNCTIONAL = "wb97x"
CHARGE = 1
MULTIPLICITY = 1

# omega is physically irrelevant here because the QED-DFT formalism uses a
# coherent-state transformation of the Pauli-Fierz Kohn-Sham equations (it is
# omega-independent). We nonetheless pin a single value everywhere -- the opt and
# freq scripts previously disagreed (0.06615 vs 0.0). 0.06615 is the value used
# for the downstream frequency-dependent QED-CCSD work, so we standardize on it.
OMEGA = 0.06615

# Psi4 thread count. Leave 0 to use Psi4's own default (typically all cores).
# When fanning cells out across processes for the embarrassingly-parallel run,
# set this per process (e.g. `export PSI4_THREADS=4`) so concurrent cells don't
# oversubscribe the machine: aim for concurrent_cells * PSI4_THREADS ~= cores.
PSI4_THREADS = int(os.environ.get("PSI4_THREADS", "0"))

# Counterpoise/Hessian unit-conversion constants (from the example freq script).
AMU_TO_AU = 1822.888486209
AU_TIME_S = 2.4188843265857e-17       # atomic unit of time, seconds
C_CM_PER_S = 2.99792458e10            # speed of light, cm/s
CM_TO_EV = 1.239841984e-4


# =============================================================================
# Cavity coupling vector
# =============================================================================

def lambda_vector_for(theta_deg, phi_deg, magnitude):
    """
    Build the cavity coupling vector for a given orientation and magnitude.

    The direction comes from the (theta, phi) spherical angles via the
    cqed_scf helper; we normalize it and scale to the requested magnitude. This
    is the *only* place a lambda vector is constructed in the campaign, so the
    optimization and the frequency analysis are guaranteed to use the same one.

    Returns
    -------
    list[float]
        The lambda vector, length 3.
    """
    direction = np.asarray(
        generate_field_vector_from_theta_and_phi(theta_deg, phi_deg),
        dtype=float,
    )
    norm = np.linalg.norm(direction)
    if norm == 0.0:
        raise ValueError(
            f"Field direction for (theta={theta_deg}, phi={phi_deg}) is zero."
        )
    unit = direction / norm
    return (magnitude * unit).tolist()


# =============================================================================
# Geometry <-> Psi4 string helpers
# =============================================================================

def psi4_geometry_string(symbols, coords_angstrom):
    """
    Build a Psi4 geometry string for the cavity calculation.

    The cavity breaks rotational invariance, so we keep the frame fixed:
    no_reorient and no_com are mandatory, and symmetry is c1.
    """
    lines = [f"{CHARGE} {MULTIPLICITY}"]
    for s, p in zip(symbols, coords_angstrom):
        lines.append(f"{s:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}")
    lines.append("units angstrom")
    lines.append("no_reorient")
    lines.append("no_com")
    lines.append("symmetry c1")
    return "\n".join(lines) + "\n"


def symbols_from_psi4_geometry(geometry_string):
    """Extract the atom symbols (in order) from a Psi4 geometry string."""
    mol = psi4.geometry(geometry_string)
    return [mol.symbol(i) for i in range(mol.natom())]


def coords_angstrom_from_opt_result(opt_result):
    """
    Convert a bfgs_optimize result's coordinate vector (flat, bohr) into an
    (natom, 3) array in angstrom.
    """
    coords_bohr = np.asarray(opt_result.x, dtype=float).reshape(-1, 3)
    return coords_bohr / ANGSTROM_TO_BOHR


# =============================================================================
# Calculator construction
# =============================================================================

def build_cqed_calculator(lambda_vector, omega=OMEGA):
    """Construct the CQEDCalculator used by both opt and freq for a cell."""
    # Both example scripts set these globally as well as passing them in; mirror
    # that so nothing depends on unset global Psi4 state (basis, grid, etc.).
    if PSI4_THREADS > 0:
        psi4.set_num_threads(PSI4_THREADS)
    psi4.set_options(PSI4_OPTIONS)
    return CQEDCalculator(
        lambda_vector=lambda_vector,
        psi4_options=PSI4_OPTIONS,
        omega=omega,
        density_fitting=True,
        charge=CHARGE,
        multiplicity=MULTIPLICITY,
        functional=FUNCTIONAL,
    )


# =============================================================================
# ASE calculator wrapping the full (unprojected) QED-DFT gradient
# =============================================================================

class CQED_DFT_Gradient(Calculator):
    """
    ASE Calculator adapter around CQEDCalculator.energy_and_gradient.

    Used only for the finite-difference Hessian. The gradients here are the
    FULL Cartesian gradients (translations/rotations are NOT projected during
    the finite differences); the rigid-body projection is applied afterward on
    the assembled Hessian, in project_and_diagonalize_hessian.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, cqed_calculator_instance):
        super().__init__()
        self.backend_calc = cqed_calculator_instance

    def calculate(self, atoms=None, properties=("energy", "forces"),
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)

        pos = atoms.get_positions()
        syms = atoms.get_chemical_symbols()

        mol_str = "\n".join(
            f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}"
            for s, p in zip(syms, pos)
        )
        mol_str += f"\n{CHARGE} {MULTIPLICITY}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com"

        energy_hartree, grad_hartree_bohr, _ = self.backend_calc.energy_and_gradient(mol_str)

        self.results["energy"] = energy_hartree * Hartree

        grad_array = np.asarray(grad_hartree_bohr).reshape(-1, 3)
        # ASE wants forces (= -gradient).
        self.results["forces"] = -grad_array * (Hartree / Bohr)


# =============================================================================
# Mass-weighted rigid-body projector and projected Hessian
# =============================================================================

def rigid_body_projector_mw(coords_bohr, masses_au, svd_tol=1e-10):
    """
    Build the mass-weighted projector that removes the 6 (or 5) rigid-body
    translation/rotation modes from a Cartesian Hessian.
    """
    coords_bohr = np.asarray(coords_bohr, dtype=float)
    masses_au = np.asarray(masses_au, dtype=float)

    natom = coords_bohr.shape[0]
    total_mass = np.sum(masses_au)

    com = np.sum(coords_bohr * masses_au[:, None], axis=0) / total_mass
    x = coords_bohr - com
    sqrtm = np.sqrt(masses_au)

    modes = []

    # Mass-weighted translations
    for axis in range(3):
        mode = np.zeros((natom, 3))
        mode[:, axis] = sqrtm
        modes.append(mode.reshape(-1))

    # Mass-weighted rotations
    for axis in range(3):
        axis_vec = np.eye(3)[axis]
        mode = np.cross(axis_vec[None, :], x)
        mode *= sqrtm[:, None]
        modes.append(mode.reshape(-1))

    B = np.column_stack(modes)

    U, s, _ = np.linalg.svd(B, full_matrices=False)
    keep = s > svd_tol
    Q = U[:, keep]

    Pmw = np.eye(3 * natom) - Q @ Q.T
    return Pmw, Q, s


def project_and_diagonalize_hessian(atoms, vib):
    """
    Extract the ASE Cartesian Hessian, convert to atomic units, mass-weight,
    project out translations/rotations, diagonalize, and return frequencies.

    Faithful to the example freq script, but returns everything in a dict with
    no hard-coded filenames so the caller controls output naming.
    """
    natom = len(atoms)

    # ASE Hessian from finite differences of forces, in eV / Angstrom^2.
    H_ase = vib.get_vibrations().get_hessian()
    H_ase = np.asarray(H_ase).reshape(3 * natom, 3 * natom)

    # Convert eV / Angstrom^2 -> Hartree / bohr^2.
    H_cart_au = H_ase * (Bohr**2 / Hartree)
    H_cart_au = 0.5 * (H_cart_au + H_cart_au.T)

    coords_bohr = atoms.get_positions() * ANGSTROM_TO_BOHR
    masses_au = atoms.get_masses() * AMU_TO_AU

    Pmw, Q, singular_values = rigid_body_projector_mw(coords_bohr, masses_au)

    sqrtm_vec = np.repeat(np.sqrt(masses_au), 3)
    Minv_sqrt = np.diag(1.0 / sqrtm_vec)

    H_mw = Minv_sqrt @ H_cart_au @ Minv_sqrt

    H_locked_mw = Pmw @ H_mw @ Pmw
    H_locked_mw = 0.5 * (H_locked_mw + H_locked_mw.T)

    evals, evecs = np.linalg.eigh(H_locked_mw)

    omega_au = np.sign(evals) * np.sqrt(np.abs(evals))
    freq_cm = omega_au / (2.0 * np.pi * C_CM_PER_S * AU_TIME_S)

    return {
        "H_cart_au": H_cart_au,
        "H_mw": H_mw,
        "H_locked_mw": H_locked_mw,
        "evals": evals,
        "evecs": evecs,
        "freq_cm": freq_cm,
        "projector_rank": int(Q.shape[1]),
        "singular_values": singular_values,
    }


def zero_point_energy(freq_cm, low_cut_cm=1.0):
    """
    Compute the ZPE from projected frequencies (cm^-1).

    Near-zero residual modes (|f| <= low_cut_cm) are dropped, and only real
    (positive) frequencies contribute. Imaginary modes are reported separately
    by the caller; they should NOT silently enter the ZPE.
    """
    freq_cm = np.asarray(freq_cm, dtype=float)

    nonzero = np.abs(freq_cm) > low_cut_cm
    real_positive = freq_cm[nonzero & (freq_cm > 0.0)]
    imaginary = freq_cm[nonzero & (freq_cm < 0.0)]

    zpe_ev = 0.5 * np.sum(real_positive) * CM_TO_EV
    zpe_hartree = zpe_ev / Hartree

    return {
        "zpe_ev": float(zpe_ev),
        "zpe_hartree": float(zpe_hartree),
        "n_real_modes": int(real_positive.size),
        "n_imaginary": int(imaginary.size),
        "imaginary_freqs_cm": imaginary.tolist(),
        "lowest_freq_cm": float(freq_cm[nonzero].min()) if nonzero.any() else None,
    }
