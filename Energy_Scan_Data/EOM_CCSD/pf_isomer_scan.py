#!/usr/bin/env python3
"""
Pauli-Fierz Hamiltonian (theta, phi) scan for three nitro-isomers.

Generalized from parse_eom_build_pf.ipynb.  For each (theta, phi) on a grid,
a lambda vector of fixed magnitude is built, the Pauli-Fierz (PF) Hamiltonian
is constructed from EOM-CCSD data for each isomer (para, ortho, meta), and the
polaritonic ground state is obtained by diagonalization.

Two flavors of the PF Hamiltonian are built and diagonalized side by side for
every (theta, phi) point and every isomer:

  - "non-CS": the bilinear coupling and dipole self-energy operators are built
    from the total (electronic + nuclear) dipole, so both nuclear and
    electronic dipoles contribute to the d operator.
  - "CS": the coherent-state transformation is applied, which shifts d by its
    ground-state expectation value <d>.  As a consequence only the
    electronic part of the dipole contributes to d.  H_el and H_ph are the
    same operators in both cases, but H_blc, H_dse, the eigenvectors, and
    therefore the decomposition expectation values differ between the two.

Outputs (fixed-width text files, written under Nel_<N>_Nph_<M>/):
  - isomer_Nel_<N>_Nph_<M>_total_energies.dat          total energy per isomer, non-CS
  - isomer_Nel_<N>_Nph_<M>_total_energies_CS.dat       total energy per isomer, CS
  - isomer_Nel_<N>_Nph_<M>_energy_decomposition.dat    E_el/E_ph/E_blc/E_dse per isomer, non-CS
  - isomer_Nel_<N>_Nph_<M>_energy_decomposition_CS.dat E_el/E_ph/E_blc/E_dse per isomer, CS
  - isomer_Nel_<N>_Nph_<M>_differences.dat             per-isomer (non-CS - CS) for the
                                                        total energy and each decomposition term

Only the configuration block below should normally need editing.  In
particular, NUM_ELECTRONIC_STATES and NUM_FOCK_STATES are the two knobs for
toggling the size of the PF Hamiltonian, and they also set the output
directory/file names automatically.
"""

import os
import numpy as np
import h5py

# ============================================================================
#                              CONFIGURATION
# ============================================================================

# --- Directory containing the isomer h5 files -------------------------------
DATA_DIR = "reexternaldecisiononsubmissiontochemicalscience"

# --- Isomer label -> h5 filename (column order in the output follows this) --
ISOMERS = {
    "Para":  "para.h5",
    "Ortho": "ortho.h5",
    "Meta":  "meta.h5",
}

# --- PF Hamiltonian size knobs ----------------------------------------------
NUM_ELECTRONIC_STATES = 49   # number of electronic states used to build H_PF
NUM_FOCK_STATES       = 10    # number of photonic Fock states

# --- Photon / coupling parameters -------------------------------------------
OMEGA            = 0.066148  # photon frequency (Hartree)
LAMBDA_MAGNITUDE = 0.1       # magnitude of the coupling (lambda) vector

# --- Nuclear dipole moment (atomic units) -----------------------------------
NUCLEAR_DIPOLE_ORTHO = np.array([-6.2324621, 14.0425864, 13.8829384])  # (x, y, z) components
NUCLEAR_DIPOLE_META = np.array([9.3777113, 13.8918947,  13.0957213])  # (x, y, z) components
NUCLEAR_DIPOLE_PARA = np.array([16.6636935, 0.2903430, 13.5397750])  # (x, y, z) components
# --- Angle grid (degrees) ---------------------------------------------------
THETA_LIST = np.linspace(0.0, 180.0, 90)
PHI_LIST   = np.linspace(0.0, 360.0, 90)

# --- Output directory / file names, derived from the size knobs above ------
# Both the CS and non-CS Hamiltonians are always built and diagonalized, so
# every run produces all five files below.
_OUT_TAG  = f"Nel_{NUM_ELECTRONIC_STATES}_Nph_{NUM_FOCK_STATES}"
_BASENAME = f"isomer_{_OUT_TAG}"
OUTPUT_DIR = _OUT_TAG

TOTAL_ENERGY_FILE    = os.path.join(OUTPUT_DIR, f"{_BASENAME}_total_energies.dat")
TOTAL_ENERGY_FILE_CS = os.path.join(OUTPUT_DIR, f"{_BASENAME}_total_energies_CS.dat")
DECOMP_FILE          = os.path.join(OUTPUT_DIR, f"{_BASENAME}_energy_decomposition.dat")
DECOMP_FILE_CS       = os.path.join(OUTPUT_DIR, f"{_BASENAME}_energy_decomposition_CS.dat")
DIFFERENCE_FILE      = os.path.join(OUTPUT_DIR, f"{_BASENAME}_differences.dat")


# ============================================================================
#                              CORE ROUTINES
# ============================================================================

def generate_lambda_vec_from_theta_and_phi(theta, phi):
    """Unit vector from spherical angles theta, phi given in DEGREES."""
    theta = theta * np.pi / 180.0
    phi   = phi   * np.pi / 180.0

    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    return np.array([x, y, z])


def parse_cq_h5_data(h5_file_path, verbose=False, nuclear_dipole=None):
    """Parse EOM-CCSD data from a ChronusQuantum h5 file.

    Returns
    -------
    reference_energy   : (1,) array
    correlation_energy : (1,) array
    excitation_energies: (num_roots,) array
    master_dipole      : (total_states, total_states, 3) array
    """
    if verbose:
        print(f"Parsing EOMCC data from: {h5_file_path}")


    if nuclear_dipole is None:
        if "ortho" in h5_file_path:
            nuclear_dipole = NUCLEAR_DIPOLE_ORTHO
            print("  Using ortho nuclear dipole")
        elif "meta" in h5_file_path:
            nuclear_dipole = NUCLEAR_DIPOLE_META
            print("  Using meta nuclear dipole")
        elif "para" in h5_file_path:
            nuclear_dipole = NUCLEAR_DIPOLE_PARA
            print("  Using para nuclear dipole")
        else:
            raise ValueError(f"Unknown isomer in file path: {h5_file_path}")

    with h5py.File(h5_file_path, "r") as f:
        # 1. Scalars
        reference_energy = f["CC/REFERENCE_ENERGY"][()]
        correlation_energy = f["CC/CORRELATION_ENERGY"][()]

        # 2. Excitation energies
        excitation_energies = f["CC/EXCITATION_ENERGIES"][:]
        num_excited_roots = len(excitation_energies)
        total_states = num_excited_roots + 1  # + ground state (index 0)

        # 3. Dipole datasets.  Transpose flips the row-major h5 layout to CQ's
        #    intended column-major layout.
        gs_dipole = f["CC/GROUND_STATE_DIPOLE"][:].T
        g2e_dipole = f["CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE"][:].T
        e2g_dipole = f["CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE"][:].T
        e2e_dipole = f["CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE"][:].T

        # Ensure the 3 Cartesian components are the LAST dimension.
        def ensure_cartesian_last(arr):
            if arr.shape[0] == 3 and arr.shape[-1] != 3:
                return np.moveaxis(arr, 0, -1)
            return arr

        gs_dipole = ensure_cartesian_last(gs_dipole)
        g2e_dipole = ensure_cartesian_last(g2e_dipole)
        e2g_dipole = ensure_cartesian_last(e2g_dipole)
        e2e_dipole = ensure_cartesian_last(e2e_dipole)

        print("  Parsed dipole shapes:")
        print(f"    Ground state: {gs_dipole.shape}")
        print(f"    G->E transitions: {g2e_dipole.shape}")
        print(f"    E->G transitions: {e2g_dipole.shape}")
        print(f"    E->E transitions: {e2e_dipole.shape}")

        # 4. Assemble the master dipole tensor (total_states, total_states, 3)
        electronic_dipole = np.zeros((total_states, total_states, 3))
        nuclear_contribution = np.zeros_like(electronic_dipole)

        electronic_dipole[0, 0, :] = gs_dipole.ravel()
        electronic_dipole[0, 1:, :] = g2e_dipole.reshape(num_excited_roots, 3)
        electronic_dipole[1:, 0, :] = e2g_dipole.reshape(num_excited_roots, 3)
        electronic_dipole[1:, 1:, :] = e2e_dipole.reshape(
            num_excited_roots, num_excited_roots, 3
        )
        # The nuclear dipole is a c-number in electronic space (mu_nuc * identity),
        # so it contributes ONLY to the diagonal <p|mu|p>, never to transition elements.
        diag = np.arange(total_states)
        nuclear_contribution[diag, diag, :] = nuclear_dipole
        electronic_dipole *= -1.0 # ChronusQuantum uses opposite sign convention
        # you need to add the nuclear dipole to the electronic dipole to get the total dipole
        total_dipole = electronic_dipole + nuclear_contribution

        print(f"   <0|mu|0> dipole is {total_dipole[0, 0, :]}")
        print(f"   <1|mu|1> dipole is {total_dipole[1, 1, :]}")
        print(f"   <2|mu|2> dipole is {total_dipole[2, 2, :]}")
        print(f"   <3|mu|3> dipole is {total_dipole[3, 3, :]}")

    return reference_energy, correlation_energy, excitation_energies, electronic_dipole, total_dipole


def build_electronic_energies(ref_E, corr_E, excit_E):
    """Build the vector of bare electronic state energies (ground + excited)."""
    dim_el = len(excit_E) + 1
    E_el = np.zeros(dim_el)
    E_el[0] = np.real(ref_E[0] + corr_E[0])
    for i in range(1, dim_el):
        E_el[i] = np.real(ref_E[0] + corr_E[0] + excit_E[i - 1])
    return E_el


def build_ladder_operator(dim):
    """Lowering and raising operators for a truncated Fock space of size dim."""
    lower_mat = np.zeros((dim, dim))
    for i in range(dim - 1):
        lower_mat[i, i + 1] = np.sqrt(i + 1)
    raise_mat = lower_mat.conj().T
    return lower_mat, raise_mat


def build_PF_Hamiltonian(dim_ph, dim_el, omega, lambda_vec, e_el, electronic_dipoles, total_dipoles):
    """Build the Pauli-Fierz Hamiltonian and its component pieces, both with
    (CS) and without (non-CS) the coherent-state transformation.

    Parameters
    ----------
    dim_ph     : number of photonic Fock states
    dim_el     : number of electronic states
    omega      : photon frequency (Hartree)
    lambda_vec : (3,) coupling vector
    e_el       : (>=dim_el,) bare electronic energies
    electronic_dipoles : (>=dim_el, >=dim_el, 3) electronic dipole tensor
    total_dipoles      : (>=dim_el, >=dim_el, 3) total dipole tensor, electronic + nuclear

    Returns
    -------
    H_total, H_total_cs, H_el_full, H_ph_full, H_blc, H_dse, H_blc_cs, H_dse_cs
        where H_el_full = kron(H_el, I_ph) and H_ph_full = kron(I_el, H_ph) are
        common to both flavors; H_blc/H_dse are the non-CS (total-dipole)
        bilinear-coupling/dipole-self-energy operators, and H_blc_cs/H_dse_cs
        are their CS (electronic-dipole-only, shifted by <d>) counterparts.
    """
    # ladder operators and identities
    a, a_dag = build_ladder_operator(dim_ph)
    I_ph = np.eye(dim_ph)
    I_el = np.eye(dim_el)

    # bare electronic Hamiltonian (truncated to dim_el states)
    H_el = np.diag(e_el[:dim_el])

    # photon Hamiltonian
    H_ph = omega * (a_dag @ a + 0.5 * np.eye(dim_ph))

    # d[p,q] = sum_i lambda_vec[i] * dipoles[p,q,i]
    d_total = np.einsum("i,pqi->pq", lambda_vec, total_dipoles)
    d_elec = np.einsum("i,pqi->pq", lambda_vec, electronic_dipoles)

    # truncate to the requested number of electronic states
    d_total = d_total[:dim_el, :dim_el]

    # coherent-state shift: subtract <0|d|0> on the diagonal (electronic dipole only)
    d_0 = d_elec[0, 0] * I_el
    d_cs = d_elec[:dim_el, :dim_el] - d_0

    # bilinear coupling
    H_blc_cs = -np.sqrt(omega / 2.0) * np.kron(d_cs, (a + a_dag))
    H_blc = -np.sqrt(omega / 2.0) * np.kron(d_total, (a + a_dag))

    # dipole self energy
    H_dse_cs = 0.5 * np.kron(d_cs @ d_cs, I_ph)
    H_dse = 0.5 * np.kron(d_total @ d_total, I_ph)

    # full pieces and total
    H_el_full = np.kron(H_el, I_ph)
    H_ph_full = np.kron(I_el, H_ph)
    H_total_cs = H_el_full + H_ph_full + H_blc_cs + H_dse_cs
    H_total = H_el_full + H_ph_full + H_blc + H_dse

    return H_total, H_total_cs, H_el_full, H_ph_full, H_blc, H_dse, H_blc_cs, H_dse_cs


# ============================================================================
#                              OUTPUT HELPERS
# ============================================================================

def _format_header(angle_cols, value_cols):
    """Build a fixed-width header line plus separator rule."""
    parts = [f"{c:>10}" for c in angle_cols[:2]]          # theta, phi
    parts += [f"{c:>13}" for c in angle_cols[2:]]          # Ex, Ey, Ez
    parts += [f"{c:>21}" for c in value_cols]              # energies
    header = "".join(parts)
    rule = "-" * len(header)
    return header + "\n" + rule + "\n"


def _format_row(theta, phi, lam, values):
    """Build one fixed-width data row."""
    parts = [f"{theta:>10.3f}", f"{phi:>10.3f}"]
    parts += [f"{v:>13.6f}" for v in lam]
    parts += [f"{v:>21.12f}" for v in values]
    return "".join(parts) + "\n"


# ============================================================================
#                                 MAIN
# ============================================================================

def main():
    labels = list(ISOMERS.keys())

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- Load and pre-process each isomer's data ---------------------------
    isomer_data = {}
    for label, fname in ISOMERS.items():
        path = os.path.join(DATA_DIR, fname)
        ref_E, corr_E, excit_E, electronic_dipole, total_dipoles = parse_cq_h5_data(path, verbose=True)

        # symmetrize the transition-dipole matrix (bra/ket swap)
        electronic_dipoles_sym = 0.5 * (electronic_dipole + np.transpose(electronic_dipole, axes=(1, 0, 2)))
        total_dipoles_sym = 0.5 * (total_dipoles + np.transpose(total_dipoles, axes=(1, 0, 2)))

        E_el = build_electronic_energies(ref_E, corr_E, excit_E)
        isomer_data[label] = {"E_el": E_el, "electronic_dipoles": electronic_dipoles_sym, "total_dipoles": total_dipoles_sym}

        if NUM_ELECTRONIC_STATES > len(E_el):
            raise ValueError(
                f"{label}: NUM_ELECTRONIC_STATES={NUM_ELECTRONIC_STATES} "
                f"exceeds available states ({len(E_el)})."
            )

    # --- Open output files and write headers -------------------------------
    angle_cols = ["theta", "phi", "Ex", "Ey", "Ez"]

    total_header = _format_header(angle_cols, [f"{l}_E" for l in labels])

    decomp_cols = []
    for l in labels:
        decomp_cols += [f"{l}_E_el", f"{l}_E_ph", f"{l}_E_blc", f"{l}_E_dse"]
    decomp_header = _format_header(angle_cols, decomp_cols)

    diff_cols = []
    for l in labels:
        diff_cols += [f"{l}_dE_total", f"{l}_dE_el", f"{l}_dE_ph", f"{l}_dE_blc", f"{l}_dE_dse"]
    diff_header = _format_header(angle_cols, diff_cols)

    f_tot    = open(TOTAL_ENERGY_FILE, "w")
    f_tot_cs = open(TOTAL_ENERGY_FILE_CS, "w")
    f_dec    = open(DECOMP_FILE, "w")
    f_dec_cs = open(DECOMP_FILE_CS, "w")
    f_dif    = open(DIFFERENCE_FILE, "w")

    f_tot.write(total_header)
    f_tot_cs.write(total_header)
    f_dec.write(decomp_header)
    f_dec_cs.write(decomp_header)
    f_dif.write(diff_header)

    # --- Scan over the angle grid ------------------------------------------
    n_points = len(THETA_LIST) * len(PHI_LIST)
    print(f"\nRunning PF scan: {n_points} (theta, phi) points x "
          f"{len(labels)} isomers x 2 Hamiltonians (non-CS, CS)")
    print(f"  electronic states = {NUM_ELECTRONIC_STATES}, "
          f"Fock states = {NUM_FOCK_STATES}, omega = {OMEGA}")
    print(f"  output directory  = {OUTPUT_DIR}")

    count = 0
    for theta in THETA_LIST:
        for phi in PHI_LIST:
            lam = generate_lambda_vec_from_theta_and_phi(theta, phi) * LAMBDA_MAGNITUDE

            totals    = {}
            totals_cs = {}
            decomp    = {}
            decomp_cs = {}
            for label in labels:
                d = isomer_data[label]
                H_total, H_total_cs, H_el_full, H_ph_full, H_blc, H_dse, H_blc_cs, H_dse_cs = build_PF_Hamiltonian(
                    dim_ph=NUM_FOCK_STATES,
                    dim_el=NUM_ELECTRONIC_STATES,
                    omega=OMEGA,
                    lambda_vec=lam,
                    e_el=d["E_el"],
                    electronic_dipoles=d["electronic_dipoles"],
                    total_dipoles=d["total_dipoles"]
                )

                # non-CS: diagonalize H_total (total-dipole H_blc/H_dse)
                eigs, C = np.linalg.eigh(H_total)
                psi0 = C[:, 0]  # polaritonic ground state, non-CS

                # CS: diagonalize H_total_cs (electronic-dipole-only, shifted H_blc_cs/H_dse_cs)
                eigs_cs, C_cs = np.linalg.eigh(H_total_cs)
                psi0_cs = C_cs[:, 0]  # polaritonic ground state, CS

                totals[label]    = eigs[0]
                totals_cs[label] = eigs_cs[0]

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

            # total energies (non-CS, CS)
            f_tot.write(_format_row(theta, phi, lam, [totals[l] for l in labels]))
            f_tot_cs.write(_format_row(theta, phi, lam, [totals_cs[l] for l in labels]))

            # decomposition (non-CS, CS), all isomers side by side
            decomp_vals = []
            decomp_cs_vals = []
            for l in labels:
                decomp_vals.extend(decomp[l])
                decomp_cs_vals.extend(decomp_cs[l])
            f_dec.write(_format_row(theta, phi, lam, decomp_vals))
            f_dec_cs.write(_format_row(theta, phi, lam, decomp_cs_vals))

            # model differences per isomer: non-CS minus CS, for total + each component
            diff_vals = []
            for l in labels:
                d_total_diff = totals[l] - totals_cs[l]
                d_el_diff, d_ph_diff, d_blc_diff, d_dse_diff = (
                    decomp[l][i] - decomp_cs[l][i] for i in range(4)
                )
                diff_vals += [d_total_diff, d_el_diff, d_ph_diff, d_blc_diff, d_dse_diff]
            f_dif.write(_format_row(theta, phi, lam, diff_vals))

            count += 1
            if count % 500 == 0 or count == n_points:
                print(f"  ... {count}/{n_points} points done")

    f_tot.close()
    f_tot_cs.close()
    f_dec.close()
    f_dec_cs.close()
    f_dif.close()

    print("\nDone. Wrote:")
    print(f"  {TOTAL_ENERGY_FILE}     (total energies, non-CS: {', '.join(labels)})")
    print(f"  {TOTAL_ENERGY_FILE_CS}  (total energies, CS: {', '.join(labels)})")
    print(f"  {DECOMP_FILE}           (decomposition, non-CS: E_el, E_ph, E_blc, E_dse per isomer)")
    print(f"  {DECOMP_FILE_CS}        (decomposition, CS: E_el, E_ph, E_blc, E_dse per isomer)")
    print(f"  {DIFFERENCE_FILE}       (non-CS minus CS: total + E_el/E_ph/E_blc/E_dse per isomer)")


if __name__ == "__main__":
    main()
