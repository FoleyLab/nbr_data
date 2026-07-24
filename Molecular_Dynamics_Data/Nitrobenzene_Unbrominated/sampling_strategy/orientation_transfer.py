"""
orientation_transfer.py
=======================

Single source of truth for turning a set of Cartesian coordinates into a proper
ORTHONORMAL molecular body frame, and for transferring a field orientation
(expressed as body-frame components) onto an intermediate's own geometry to build
the cavity coupling vector lambda.

This is deliberately independent of cqed-scf so it can be unit-tested cheaply.
In the real tree it should live inside (or be imported by) cavity_common.py.

Frame convention (matches NitrobenzeneOrientation, made orthonormal & gauge-safe)
--------------------------------------------------------------------------------
  x_hat  : ipso-carbon -> N  (the C-N bond direction)
  z_hat  : ring normal from a carbon triplet cross product, orthogonalized
           against x_hat (Gram-Schmidt). Its SIGN is molecular, fixed by the
           carbon atom ORDERING, *not* by the lab z axis. NitrobenzeneOrientation
           forced z[2]>0 (a lab flip); we do NOT, because when the molecule turns
           over that flip silently relabels the two ring faces -- which matters
           the moment Br breaks the face symmetry in an intermediate.
  y_hat  : z_hat x x_hat   (completes a right-handed orthonormal frame)

A field unit vector f has body components (a,b,c) = (x_hat.f, y_hat.f, z_hat.f).
Because the frame is orthonormal this is fully invertible -- unlike the two
arccos projections (theta,phi), which throw away the sign of b.

IMPORTANT ATOM-ORDERING ASSUMPTION
----------------------------------
`ring_triplet` defaults to the first three carbons in FILE ORDER, exactly as
NitrobenzeneOrientation did. For the transfer to be meaningful, every
intermediate .xyz must list its atoms so that (a) the ipso carbon is identified
the same way and (b) the first three carbons define the ring normal with the
same molecular handedness as the reactant. build_intermediate_frame() asserts
Br-on-minus-z and near-orthogonality so a mismatch fails LOUDLY instead of
producing a quietly wrong lambda.
"""

import numpy as np


def _unit(v):
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("degenerate vector in frame construction")
    return v / n


def find_ipso_carbon(symbols, coords):
    """Carbon nearest the (single) nitrogen -- the ipso carbon."""
    C = [i for i, s in enumerate(symbols) if s == "C"]
    N = [i for i, s in enumerate(symbols) if s == "N"][0]
    d = np.linalg.norm(coords[C] - coords[N], axis=1)
    return C[int(np.argmin(d))], N


def body_frame(symbols, coords, ring_triplet=None):
    """
    Return (x_hat, y_hat, z_hat): an orthonormal molecular frame.

    ring_triplet : optional 3 carbon indices used for the ring normal.
                   Defaults to the first three carbons in file order
                   (matching NitrobenzeneOrientation._choose_ring_triplet).
    """
    coords = np.asarray(coords, float)
    ipso, N = find_ipso_carbon(symbols, coords)
    if ring_triplet is None:
        C = [i for i, s in enumerate(symbols) if s == "C"]
        ring_triplet = C[:3]
    i, j, k = ring_triplet

    x_hat = _unit(coords[N] - coords[ipso])                 # ipso -> N
    z_raw = np.cross(coords[j] - coords[i], coords[k] - coords[i])
    z_hat = _unit(z_raw - (z_raw @ x_hat) * x_hat)          # orthogonalize vs x
    y_hat = np.cross(z_hat, x_hat)                          # right-handed
    return x_hat, y_hat, z_hat


def field_body_components(symbols, coords, field_vector, ring_triplet=None):
    """(a,b,c) = components of the (normalized) field in the body frame."""
    x_hat, y_hat, z_hat = body_frame(symbols, coords, ring_triplet)
    f = _unit(np.asarray(field_vector, float))
    return np.array([x_hat @ f, y_hat @ f, z_hat @ f])


def build_intermediate_frame(symbols, coords, ring_triplet=None,
                             br_symbol="Br", check=True):
    """
    Body frame for an intermediate, with safety checks that catch the silent
    frame-mismatch failure mode.

    Asserts (when check=True):
      * the frame is orthonormal to ~1e-6
      * Br sits on the -z side (the reference convention). If Br is on +z the
        atom ordering / handedness disagrees with the reactant and the transfer
        would be wrong -- we raise rather than proceed.
    Returns (x_hat, y_hat, z_hat).
    """
    coords = np.asarray(coords, float)
    x_hat, y_hat, z_hat = body_frame(symbols, coords, ring_triplet)

    if check:
        M = np.vstack([x_hat, y_hat, z_hat])
        ortho_err = np.abs(M @ M.T - np.eye(3)).max()
        if ortho_err > 1e-6:
            raise AssertionError(f"frame not orthonormal (err={ortho_err:.1e})")

        br = [idx for idx, s in enumerate(symbols) if s == br_symbol]
        if br:
            ring_C = [i for i, s in enumerate(symbols) if s == "C"]
            centroid = coords[ring_C].mean(0)
            c_br = (coords[br[0]] - centroid) @ z_hat
            if c_br > 0:
                raise AssertionError(
                    f"Br is on +z (z.(Br-centroid)={c_br:+.3f}) but convention is "
                    f"Br on -z. Atom ordering / ring-triplet handedness disagrees "
                    f"with the reactant frame -- fix before running the campaign.")
    return x_hat, y_hat, z_hat


def lambda_for_intermediate(abc, symbols, coords, magnitude,
                            ring_triplet=None, check=True):
    """
    Build the Cartesian cavity coupling vector lambda for an intermediate whose
    geometry is `coords`, given a field orientation `abc=(a,b,c)` expressed in
    the (reactant) body frame.

        lambda = magnitude * (a*x_hat + b*y_hat + c*z_hat)

    where x_hat,y_hat,z_hat are the intermediate's OWN body axes. Because both
    reactant and intermediate frames are built identically, the body components
    (a,b,c) denote the same physical orientation of the field relative to the
    molecule in each. This is the generalization of lambda_vector_for(theta,phi)
    to a full, sign-complete body vector.
    """
    a, b, c = abc
    x_hat, y_hat, z_hat = build_intermediate_frame(
        symbols, coords, ring_triplet=ring_triplet, check=check)
    direction = a * x_hat + b * y_hat + c * z_hat
    # |abc| should be ~1; renormalize defensively so |lambda| == magnitude exactly
    return magnitude * direction / np.linalg.norm(direction)


if __name__ == "__main__":
    # tiny self-test on the reactant frame 0: recover the lab field from (a,b,c)
    import json, re
    path = "/mnt/user-data/uploads/nitrobenzene_direction_A_wb97x_d_4000_ts.xyz"
    L = open(path).read().splitlines()
    nat = int(L[0]); syms = [L[2 + j].split()[0] for j in range(nat)]
    xyz = np.array([[float(x) for x in L[2 + j].split()[1:4]] for j in range(nat)])
    f_hat = np.array([0.78781234, 0.55163218, 0.27395924])
    abc = field_body_components(syms, xyz, f_hat)
    x_hat, y_hat, z_hat = body_frame(syms, xyz)
    f_rec = abc[0] * x_hat + abc[1] * y_hat + abc[2] * z_hat
    print("abc =", abc, " |abc| =", np.linalg.norm(abc))
    print("field recovery error:", np.linalg.norm(f_rec - f_hat))
