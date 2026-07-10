#!/usr/bin/env python
"""
Excited-state permanent dipoles and (all) transition dipoles for the low-lying
singlet states of the ortho bromonitro-arenium (Wheland) cation.

Method:  ADC(2) on top of an RHF reference, driven through the `adcc` package
         using a Psi4 SCF wavefunction.

Why ADC(2) rather than plain LR-TDDFT/CIS?
  You asked for THREE kinds of quantities:
     (1) permanent dipole of each state (S0, S1, S2, ...)
     (2) ground -> excited transition dipoles (0 -> n)
     (3) excited -> excited transition dipoles (n -> m, "state-to-state")
  Native Psi4 TD-DFT (`tdscf_excitations`) only gives you (2) cleanly.
  adcc gives all three out of the box:
     - state.state_dipole_moment       -> (1) for the excited states
     - state.transition_dipole_moment  -> (2)
     - adcc.State2States(...)          -> (3)
  ADC(2) is still a "cheap correlated" method. If you want something even
  cheaper and essentially CIS-like, use adcc.adc1(...) (ADC(1) == CIS).
  A pure LR-TDDFT alternative is sketched at the bottom of this file.

Requirements:
  conda install psi4 -c conda-forge
  conda install adcc -c conda-forge      # or: pip install adcc

Run with:
  python ortho_excited_dipoles.py
"""

import numpy as np
import psi4
import adcc

# ----------------------------------------------------------------------------
# 0. Housekeeping
# ----------------------------------------------------------------------------
AU2DEBYE = 2.5417464519          # 1 a.u. of dipole moment = 2.5417... Debye
AU2EV    = 27.211386245988

psi4.set_memory("16 GB")         # bump this up if the ADC step runs out of memory
psi4.set_num_threads(4)
adcc.set_n_threads(4)
psi4.core.set_output_file("ortho.out", False)

# Number of SINGLET excited states to solve for.
# You asked for the "first 3 electronic states". Below we solve for 3 EXCITED
# singlets (S1, S2, S3) and also report the ground state (S0). That way you
# have both readings covered: {S0,S1,S2} or {S1,S2,S3}. Change if you like.
N_EXC = 3

# ----------------------------------------------------------------------------
# 1. Geometry  (verbatim from your input)
#    Note: "1 1" => charge = +1, multiplicity = 1 (closed-shell singlet).
#    This is a sigma-complex / Wheland cation (the sp3 C bears both an H and
#    the Br out of plane), so an RHF reference is appropriate.
#    Br (Z=35) is all-electron in the def2 basis family (no ECP needed).
#    If the "BR" label ever fails to parse, change it to "Br".
# ----------------------------------------------------------------------------
mol = psi4.geometry("""
C           -1.804928163307     1.957993763262     0.703312273806
C           -0.379708783307     1.994122833262     0.698532703806
C            0.296125016693     0.817793533262     0.710271493806
C           -2.520286433307     0.755089873262     0.736288843806
H           -2.344947113307     2.899196893262     0.691895063806
H            0.158564066693     2.933869823262     0.699142733806
H           -3.601954283307     0.764862203262     0.746931053806
N            1.767881836693     0.820900013262     0.771891313806
O            2.315054046693    -0.296733496738     0.879853723806
O            2.340645916693     1.923356243262     0.711986073806
C           -1.829967733307    -0.442167236738     0.756258983806
H           -2.356763623307    -1.389967436738     0.789740873806
C           -0.361341153307    -0.491572936738     0.714148383806
H            0.119338216693    -1.238105076738     1.350400383806
BR          -0.151212663307    -1.224162306738    -1.170925976194
1 1
units angstrom
no_reorient
no_com
symmetry c1
""")

# ----------------------------------------------------------------------------
# 2. SCF reference
# ----------------------------------------------------------------------------
psi4.set_options({
    "basis":         "6-31G",   # try def2-svpd / def2-tzvp for better excited states
    "reference":     "rhf",
    "scf_type":      "df",
    "e_convergence": 1e-8,
    "d_convergence": 1e-8,
    "maxiter":       200,
})

scf_e, wfn = psi4.energy("scf", return_wfn=True)

# S0 permanent dipole from the SCF density (atomic units, length gauge).
mu_s0 = np.asarray(psi4.variable("SCF DIPOLE"))     # shape (3,)

# ----------------------------------------------------------------------------
# 3. ADC(2) excited states
# ----------------------------------------------------------------------------
# frozen_core=True freezes the 1s (etc.) cores => cheaper, negligible effect on
# valence excitations. Drop it if you want everything correlated.
#state = adcc.adc2(wfn, n_singlets=N_EXC, frozen_core=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]) #, frozen_core=True)
# 3. ADC(2) excited states
# ----------------------------------------------------------------------------
state = adcc.adc2(
    wfn, 
    n_singlets=N_EXC, 
    frozen_core=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],
    max_iter=150,        # <--- Increase the maximum number of iterations here
)

print(state.describe())

exc_en_ev = state.excitation_energy * AU2EV          # (N_EXC,)
osc       = state.oscillator_strength                # (N_EXC,)
tdm_0n    = np.atleast_2d(state.transition_dipole_moment)   # (N_EXC,3) 0->n, a.u.
mu_exc    = np.atleast_2d(state.state_dipole_moment)        # (N_EXC,3) permanent, a.u.

# ----------------------------------------------------------------------------
# 4. Excited -> excited (state-to-state) transition dipoles
#    State2States(state, initial=i) gives transitions from excited state i to
#    every excited state j > i.  Excited index 0 == S1, index 1 == S2, ...
# ----------------------------------------------------------------------------
s2s_tdm = {}   # s2s_tdm[i] = array over final states j>i, each a length-3 dipole (a.u.)
for i in range(N_EXC - 1):
    sts = adcc.State2States(state, initial=i)
    s2s_tdm[i] = np.atleast_2d(sts.transition_dipole_moment)

# ----------------------------------------------------------------------------
# 5. Collect everything into simple state-indexed containers.
#    Global state index: 0 = S0 (ground), 1..N_EXC = S1..S(N_EXC).
# ----------------------------------------------------------------------------
n_states = N_EXC + 1
labels   = ["S%d" % k for k in range(n_states)]

# Permanent dipoles per state (a.u.)
perm_dipole = np.zeros((n_states, 3))
perm_dipole[0] = mu_s0
for n in range(N_EXC):
    perm_dipole[n + 1] = mu_exc[n]

# Full transition-dipole "matrix" (a.u.), symmetric in magnitude.
tdm = np.zeros((n_states, n_states, 3))
# ground -> excited
for n in range(N_EXC):
    tdm[0, n + 1] = tdm_0n[n]
    tdm[n + 1, 0] = tdm_0n[n]
# excited -> excited
for i in range(N_EXC - 1):
    for k, j in enumerate(range(i + 1, N_EXC)):
        d = s2s_tdm[i][k]
        tdm[i + 1, j + 1] = d
        tdm[j + 1, i + 1] = d

# ----------------------------------------------------------------------------
# 6. Report (printed to stdout and written to a text file)
# ----------------------------------------------------------------------------
def fmt_vec(v):
    return "[%9.4f %9.4f %9.4f]" % (v[0], v[1], v[2])

lines = []
lines.append("=" * 78)
lines.append("ADC(2)/def2-SVP  --  ortho bromonitro-arenium cation (charge +1, singlet)")
lines.append("=" * 78)
lines.append("")
lines.append("Vertical excitation energies")
lines.append("-" * 78)
lines.append("  state      E_exc (eV)     osc. strength")
for n in range(N_EXC):
    lines.append("   S%-3d      %10.4f       %10.5f" % (n + 1, exc_en_ev[n], osc[n]))
lines.append("")

lines.append("Permanent dipole moments  (state |mu|)")
lines.append("-" * 78)
lines.append("  state        mu (a.u. vector)                |mu| (a.u.)")
for s in range(n_states):
    mu_d = perm_dipole[s] #* AU2DEBYE
    lines.append("  %-4s   %-34s   %8.4f"
                 % (labels[s], fmt_vec(mu_d), np.linalg.norm(mu_d)))
lines.append("")

lines.append("Transition dipole moments  (all state pairs)")
lines.append("-" * 78)
lines.append("  pair         mu_trans (a.u. vector)          |mu| (a.u.)")
for a in range(n_states):
    for b in range(a + 1, n_states):
        d_au = tdm[a, b]
        d_d  = d_au #* AU2DEBYE
        lines.append("  %s-%s   %-34s   %8.4f"
                     % (labels[a], labels[b], fmt_vec(d_d), np.linalg.norm(d_d)))
lines.append("")
lines.append("Notes:")
lines.append("  * All vectors are length-gauge, in the input (unrotated) frame")
lines.append("    (no_com / no_reorient are respected).")
lines.append("  * S0 permanent dipole is the SCF value; excited-state permanent")
lines.append("    dipoles and all transition dipoles are ADC(2) unrelaxed values.")
lines.append("  * 1 a.u. dipole = %.6f Debye." % AU2DEBYE)

report = "\n".join(lines)
print(report)
with open("ortho_dipoles_summary.txt", "w") as fh:
    fh.write(report + "\n")

# ============================================================================
# APPENDIX: pure LR-TDDFT alternative (ground->excited transition dipoles only)
# ----------------------------------------------------------------------------
# Native Psi4 does NOT give excited-state permanent dipoles or state-to-state
# transition dipoles, so it does not fully answer your request -- but if you
# only need 0->n transition dipoles / oscillator strengths at a DFT level:
#
#   psi4.set_options({"basis": "def2-svp", "reference": "rhf",
#                     "scf_type": "df", "tdscf_states": 3, "tdscf_tda": False})
#   e, wfn = psi4.energy("td-b3lyp", return_wfn=True)      # or "td-wb97x-d3"
#   res = psi4.procrouting.response.scf_response.tdscf_excitations(
#             wfn, states=3)
#   for r in res:
#       print(r["EXCITATION ENERGY"],
#             r["ELECTRIC DIPOLE TRANSITION MOMENT (LEN)"],
#             r["OSCILLATOR STRENGTH (LEN)"])
#
# Setting tdscf_tda=True with an HF reference reproduces CIS.
# ============================================================================
