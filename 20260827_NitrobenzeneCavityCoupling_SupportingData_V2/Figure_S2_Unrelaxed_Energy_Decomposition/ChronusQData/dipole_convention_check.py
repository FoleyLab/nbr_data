#!/usr/bin/env python3
"""
Numerically test which dipole-moment convention ortho.h5 (ChronusQ EOM-CCSD)
uses, by comparing against the SCF (ground state) and ADC2 (excited/transition
state) reference data in ortho_dipole_data.txt.

Hypotheses tested for the ground/excited *state* dipoles (transition dipoles
are always purely electronic in both codes, so they are compared directly):

  H1: h5_state_dipole  ==  Total   (electronic + nuclear), same sign convention
  H2: h5_state_dipole  ==  Electronic only, same sign convention
  H3: h5_state_dipole  == -Electronic only (sign-flipped electronic)
  H4: h5_state_dipole  ==  Nuclear only
  H5: h5_state_dipole - Nuclear  ==  Electronic  (i.e. h5 = Total, but using a
      possibly different nuclear-dipole reference than the one hardcoded in
      pf_isomer_scan.py)

For each hypothesis we report the vector residual, its norm relative to the
reference-vector norm, and the cosine similarity (to separate "wrong
magnitude" from "wrong direction").
"""

import h5py
import numpy as np

DATA_DIR = "." #"reexternaldecisiononsubmissiontochemicalscience"

# ---- SCF reference (ground state), from ortho_dipole_data.txt --------------
SCF_ELECTRONIC = np.array([2.4068371, -12.9048307, -12.6523611])
SCF_NUCLEAR = np.array([-6.2324621, 14.0425864, 13.8829384])
SCF_TOTAL = SCF_ELECTRONIC + SCF_NUCLEAR  # == [-3.8256..., 1.1378..., 1.2306...]

# ---- ADC2 reference (state dipoles already fold nuclear + electronic) ------
ADC2_STATE = {
    "S0": np.array([-3.8256, 1.1378, 1.2306]),
    "S1": np.array([-0.5446, 0.5907, 0.9041]),
    "S2": np.array([-0.7950, 1.0304, 0.7518]),
    "S3": np.array([0.1065, 1.1424, 1.3641]),
}

# ---- ADC2 transition dipoles (purely electronic in both codes) ------------
ADC2_TRANS = {
    ("S0", "S1"): np.array([-0.0518, -0.2609, -0.1196]),
    ("S0", "S2"): np.array([-0.1303, -0.2721, -0.1379]),
    ("S0", "S3"): np.array([0.1213, 0.2015, 0.0903]),
    ("S1", "S2"): np.array([-1.2170, -0.8745, -0.8119]),
    ("S1", "S3"): np.array([0.7749, -0.0651, 0.5320]),
    ("S2", "S3"): np.array([0.8287, 0.9386, 0.4993]),
}


def cos_sim(a, b):
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))


def rel_resid(a, b):
    """||a-b|| / ||b||"""
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def report_state(label, h5_vec):
    print(f"\n=== {label} state dipole (h5 raw vector: {h5_vec}, |mu|={np.linalg.norm(h5_vec):.4f}) ===")
    hyps = {
        "H1 Total (same sign)":            SCF_TOTAL if label == "S0" else ADC2_STATE[label],
        "H2 Electronic (same sign)":       SCF_ELECTRONIC if label == "S0" else (ADC2_STATE[label] - SCF_NUCLEAR),
        "H3 -Electronic (sign-flipped)":  -SCF_ELECTRONIC if label == "S0" else -(ADC2_STATE[label] - SCF_NUCLEAR),
        "H4 Nuclear only":                 SCF_NUCLEAR,
    }
    for name, ref in hyps.items():
        rr = rel_resid(h5_vec, ref)
        cs = cos_sim(h5_vec, ref)
        print(f"  {name:28s} ref={np.array2string(ref, precision=4):32s} "
              f"|ref|={np.linalg.norm(ref):7.4f}  rel.resid={rr:6.3f}  cos_sim={cs:6.3f}")
    # H5: subtract nuclear from h5 and compare to electronic
    h5_minus_nuc = h5_vec - SCF_NUCLEAR
    ref_elec = SCF_ELECTRONIC if label == "S0" else (ADC2_STATE[label] - SCF_NUCLEAR)
    rr = rel_resid(h5_minus_nuc, ref_elec)
    cs = cos_sim(h5_minus_nuc, ref_elec)
    print(f"  H5 (h5 - Nuclear) vs Electronic: resid={np.array2string(h5_minus_nuc - ref_elec, precision=4)} "
          f"rel.resid={rr:6.3f} cos_sim={cs:6.3f}")


def report_trans(label_pair, h5_vec):
    ref = ADC2_TRANS[label_pair]
    mag_h5 = np.linalg.norm(h5_vec)
    mag_ref = np.linalg.norm(ref)
    cs = cos_sim(h5_vec, ref)
    print(f"  {label_pair[0]}-{label_pair[1]:3s}  h5={np.array2string(h5_vec, precision=4):26s} |mu|_h5={mag_h5:7.4f}  "
          f"ADC2={np.array2string(ref, precision=4):26s} |mu|_ADC2={mag_ref:7.4f}  "
          f"|mu| ratio={mag_h5/mag_ref:6.3f}  cos_sim={cs:6.3f}")


def main():
    path = f"{DATA_DIR}/ortho.h5"
    with h5py.File(path, "r") as f:
        gs = f["CC/GROUND_STATE_DIPOLE"][:]
        g2e = f["CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE"][:]
        e2g = f["CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE"][:]
        e2e = f["CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE"][:]  # (3, nroots, nroots)

    print("############################################################")
    print("STATE (permanent) DIPOLES")
    print("############################################################")
    report_state("S0", gs)
    report_state("S1", e2e[:, 0, 0])
    report_state("S2", e2e[:, 1, 1])
    report_state("S3", e2e[:, 2, 2])

    print("\n############################################################")
    print("TRANSITION DIPOLES  (phase of each excited state is arbitrary,")
    print("so only |mu| and NOT sign/direction is meaningful for comparison)")
    print("############################################################")
    print("\n-- symmetrized bra/ket average (as pf_isomer_scan.py does) --")
    g2e_sym = 0.5 * (g2e[:, 0:3] + e2g[:, 0:3])
    report_trans(("S0", "S1"), g2e_sym[:, 0])
    report_trans(("S0", "S2"), g2e_sym[:, 1])
    report_trans(("S0", "S3"), g2e_sym[:, 2])
    report_trans(("S1", "S2"), 0.5 * (e2e[:, 0, 1] + e2e[:, 1, 0]))
    report_trans(("S1", "S3"), 0.5 * (e2e[:, 0, 2] + e2e[:, 2, 0]))
    report_trans(("S2", "S3"), 0.5 * (e2e[:, 1, 2] + e2e[:, 2, 1]))

    print("\n-- raw (unsymmetrized) g2e vs e2g, to see bra/ket asymmetry --")
    for i, pair in enumerate([("S0", "S1"), ("S0", "S2"), ("S0", "S3")]):
        print(f"  {pair}: g2e={g2e[:, i]}  e2g={e2g[:, i]}  "
              f"|g2e|={np.linalg.norm(g2e[:, i]):.4f}  |e2g|={np.linalg.norm(e2g[:, i]):.4f}")


if __name__ == "__main__":
    main()
