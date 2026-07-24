#!/usr/bin/env python3
"""
build_manifest.py
=================

Turn an annotated QED-DFT MD trajectory (nitrobenzene reactant tumbling under a
lab-fixed cavity field) into a manifest of ~K design ORIENTATIONS at which to run
the optimize -> gate -> Hessian -> ZPE pipeline for each Wheland intermediate
(ortho, meta, para).

Design philosophy (see chat handoff)
------------------------------------
The expensive pipeline never touches MD frames; it computes the three
intermediates at chosen field orientations. So we:

  STAGE 1 -- PLACEMENT (dwell-BLIND).  The trajectory only defines the *support*
      (which orientations the molecule visits). We place design points for
      COVERAGE of that support (farthest-point sampling on great-circle
      distance), and explicitly PIN one design point at each named PES feature
      (basins + hills). This is deliberately not dwell-proportional: because this
      is a zero-rotational-KE microcanonical run, dwell time piles up at the
      high-E turning points (caustic) and starves the low-E basins, so
      dwell-weighted placement would undersample exactly the regions we care
      about. Coverage placement cannot miss a region.

  STAGE 2 -- WEIGHTING (dwell fully restored, done later in analyze_surface.py).
      Every frame's dwell weight is preserved here so the fitted surface can be
      reweighted by the true occupancy afterwards -- under microcanonical dwell,
      a Boltzmann weight, or flat -- without ever having let dwell bias placement.

Symmetry fold
-------------
E(lambda)=E(-lambda) (global parity) maps visited points to the opposite ring
face, which this trajectory never visits -> no-op here, machinery kept general.
The site mirror sigma_v (the plane containing C-N and the ring normal) maps
b -> -b and swaps the two ortho sites (and the two meta sites), leaving para
fixed. Nitrobenzene is symmetric under it, so orientations (a,b,c) and (a,-b,c)
are equivalent. We therefore fold to b>=0: one computed intermediate at
(a,|b|,c) represents BOTH mirror orientations, ~halving the QC.

Outputs
-------
  manifest.json        design orientations (body-frame vector + angles + feature
                       label + dwell weight + provenance) and campaign contract.
  frames.npz           per-frame (abc, E, dwell, path coord s, assignment) for
                       Stage-2 reweighting.
  manifest_check.png   verification: design points over the dwell map and E(s).

Usage
-----
  python build_manifest.py TRAJ.xyz [-k 25] [--out-dir .] [--seed 0]
      [--lambda-mag 0.1] [--intermediates ortho meta para]
"""

import argparse, json, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from orientation_transfer import body_frame

TIMESTEP_AU = 25.0            # atomic units of time between frames
AU_TIME_FS = 0.02418884      # 1 a.u. of time in fs


# ----------------------------------------------------------------------------- IO
def parse_trajectory(path):
    lines = open(path).read().splitlines()
    pat = re.compile(r"Step\s+(\d+)\s+E=([-\d.]+)\s+phi=([-\d.]+)\s+theta=([-\d.]+)")
    steps, E, phi, theta, coords = [], [], [], [], []
    i, syms = 0, None
    while i < len(lines):
        nat = int(lines[i].strip())
        m = pat.search(lines[i + 1])
        steps.append(int(m.group(1))); E.append(float(m.group(2)))
        phi.append(float(m.group(3))); theta.append(float(m.group(4)))
        s = [lines[i + 2 + j].split() for j in range(nat)]
        if syms is None:
            syms = [r[0] for r in s]
        coords.append([[float(x) for x in r[1:4]] for r in s])
        i += 2 + nat
    return (syms, np.array(steps), np.array(E), np.array(phi),
            np.array(theta), np.array(coords))


# ----------------------------------------------------------------- geometry / fold
def field_components_all(symbols, coords, f_hat):
    """(a,b,c) per frame in the orthonormal body frame."""
    n = len(coords)
    abc = np.empty((n, 3))
    for t in range(n):
        x_hat, y_hat, z_hat = body_frame(symbols, coords[t])
        abc[t] = [x_hat @ f_hat, y_hat @ f_hat, z_hat @ f_hat]
    return abc


def fit_lab_field(symbols, coords, phi_a, theta_a):
    """
    Recover the fixed lab-frame field direction by least squares against the
    annotated angles (validates the annotations at the same time). Uses the RAW
    (non-orthogonalized) axes that NitrobenzeneOrientation used for the angles.
    """
    C = [i for i, s in enumerate(symbols) if s == "C"]
    N = [i for i, s in enumerate(symbols) if s == "N"][0]
    d0 = np.linalg.norm(coords[0][C] - coords[0][N], axis=1)
    ipso = C[int(np.argmin(d0))]; ring = C[:3]
    n = len(coords)
    xh = np.empty((n, 3)); zh = np.empty((n, 3))
    for t in range(n):
        xv = coords[t][N] - coords[t][ipso]; xh[t] = xv / np.linalg.norm(xv)
        i, j, k = ring
        zv = np.cross(coords[t][j] - coords[t][i], coords[t][k] - coords[t][i])
        zv /= np.linalg.norm(zv)
        if zv[2] < 0:
            zv = -zv          # replicate the lab flip used for the annotations
        zh[t] = zv
    A = np.vstack([xh, zh])
    b = np.concatenate([np.cos(np.radians(phi_a)), np.cos(np.radians(theta_a))])
    f, *_ = np.linalg.lstsq(A, b, rcond=None)
    f_hat = f / np.linalg.norm(f)
    phi_pred = np.degrees(np.arccos(np.clip(xh @ f_hat, -1, 1)))
    theta_pred = np.degrees(np.arccos(np.clip(zh @ f_hat, -1, 1)))
    resid = max(np.abs(phi_pred - phi_a).max(), np.abs(theta_pred - theta_a).max())
    return f_hat, resid


def fold_site_mirror(abc):
    """Fold the sigma_v site mirror: b -> |b| (b>=0)."""
    folded = abc.copy()
    folded[:, 1] = np.abs(folded[:, 1])
    return folded


# ------------------------------------------------------------- great-circle metric
def gc_dist(u, V):
    """Great-circle distance (rad) from unit vector u to each row of unit V."""
    return np.arccos(np.clip(V @ u, -1.0, 1.0))


# --------------------------------------------------------------- feature detection
def path_coordinate(folded):
    """PC1 of the folded cloud -- the dominant swing coordinate s."""
    c0 = folded - folded.mean(0)
    _, S, Vt = np.linalg.svd(c0, full_matrices=False)
    s = c0 @ Vt[0]
    var = (S ** 2 / np.sum(S ** 2))
    return s, Vt[0], var


def detect_features(s, E, nbins=60, smooth=3):
    """
    Locate the PES features along s: the 2 deepest basins (minima of E) and the
    3 highest hills (maxima), labelled major/major/minor by height. Returns a
    list of dicts with the representative FRAME index at each feature.
    """
    edges = np.linspace(s.min(), s.max(), nbins + 1)
    cent = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(s, edges) - 1, 0, nbins - 1)
    Eb = np.array([E[idx == k].mean() if (idx == k).any() else np.nan
                   for k in range(nbins)])
    # fill gaps then light smoothing
    good = ~np.isnan(Eb)
    Eb = np.interp(np.arange(nbins), np.where(good)[0], Eb[good])
    ker = np.ones(smooth) / smooth
    Es = np.convolve(Eb, ker, mode="same")

    def local(sign):
        out = []
        for k in range(1, nbins - 1):
            if sign * Es[k] > sign * Es[k - 1] and sign * Es[k] > sign * Es[k + 1]:
                out.append(k)
        return out

    minima = sorted(local(-1), key=lambda k: Es[k])          # lowest first
    maxima = sorted(local(+1), key=lambda k: -Es[k])         # highest first

    feats = []
    for k in minima[:2]:
        feats.append(("basin", cent[k]))
    hill_labels = ["major_hill", "major_hill", "minor_hill"]
    for lab, k in zip(hill_labels, maxima[:3]):
        feats.append((lab, cent[k]))

    # representative frame: within +-1 bin-width of feature s, the frame with the
    # most extreme E (min for basin, max for hill)
    bw = cent[1] - cent[0]
    result = []
    for lab, s0 in feats:
        win = np.abs(s - s0) <= bw
        cand = np.where(win)[0]
        if len(cand) == 0:
            cand = np.array([int(np.argmin(np.abs(s - s0)))])
        pick = cand[np.argmin(E[cand])] if lab == "basin" else cand[np.argmax(E[cand])]
        result.append({"label": lab, "s": float(s0), "frame": int(pick)})
    return result, (cent, Es)


# --------------------------------------------------------------- coverage placement
def farthest_point(folded, k, seed_frames):
    """
    Farthest-point (maximin) sampling on great-circle distance, seeded with the
    pinned feature frames. Density-blind: guarantees coverage, cannot miss a
    region of the support.
    """
    selected = list(dict.fromkeys(seed_frames))            # dedupe, keep order
    dmin = np.full(len(folded), np.inf)
    for idx in selected:
        dmin = np.minimum(dmin, gc_dist(folded[idx], folded))
    while len(selected) < k:
        nxt = int(np.argmax(dmin))
        if dmin[nxt] <= 0:
            break                                          # support exhausted
        selected.append(nxt)
        dmin = np.minimum(dmin, gc_dist(folded[nxt], folded))
    return selected


def voronoi_weights(folded, design_frames):
    """Assign every frame to nearest design point; dwell weight = time fraction."""
    D = np.stack([gc_dist(folded[i], folded) for i in design_frames])  # (K, n)
    assign = np.argmin(D, axis=0)
    counts = np.bincount(assign, minlength=len(design_frames))
    return assign, counts / counts.sum()


# ------------------------------------------------------------------------ manifest
def make_id(rank, abc):
    theta = np.degrees(np.arccos(np.clip(abc[2], -1, 1)))
    phi = np.degrees(np.arctan2(abc[1], abc[0]))
    return f"ori_{rank:02d}_t{theta:04.0f}_p{phi:+04.0f}".replace("+", "p").replace("-", "m")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory")
    ap.add_argument("-k", type=int, default=25, help="number of design orientations")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--lambda-mag", type=float, default=0.1)
    ap.add_argument("--intermediates", nargs="+",
                    default=["ortho", "meta", "para"])
    args = ap.parse_args()
    np.random.seed(args.seed)

    syms, steps, E, phi_a, theta_a, coords = parse_trajectory(args.trajectory)
    n = len(steps)
    Ha2kcal = 627.5094740631
    E_kcal = (E - E.min()) * Ha2kcal

    f_hat, resid = fit_lab_field(syms, coords, phi_a, theta_a)
    print(f"[frame] {n} frames; lab field f_hat={f_hat.round(5).tolist()}; "
          f"annotation residual {resid:.4f} deg")

    abc = field_components_all(syms, coords, f_hat)
    folded = fold_site_mirror(abc)

    s, pc1, var = path_coordinate(folded)
    if np.corrcoef(s, E)[0, 1] < 0:
        s, pc1 = -s, -pc1
    print(f"[path] PC1 variance fraction {var[0]:.2f}")

    features, (cent, Es) = detect_features(s, E_kcal)
    feat_frames = [f["frame"] for f in features]
    print("[features] " + ", ".join(f"{f['label']}@s={f['s']:+.2f}(frame {f['frame']})"
                                     for f in features))

    design_frames = farthest_point(folded, args.k, feat_frames)
    assign, weights = voronoi_weights(folded, design_frames)
    print(f"[placement] {len(design_frames)} design orientations "
          f"({len(feat_frames)} feature-pinned + coverage fill)")

    # per-frame dwell (uniform in time), pack for Stage 2
    dwell = np.full(n, 1.0 / n)
    np.savez(f"{args.out_dir}/frames.npz",
             abc=abc, folded=folded, E_hartree=E, s=s, dwell=dwell,
             assignment=assign, design_frames=np.array(design_frames),
             pc1=pc1, f_hat=f_hat)

    # order design points along s for readable IDs; feature ones keep their label
    order = np.argsort([s[i] for i in design_frames])
    feat_of = {f["frame"]: f["label"] for f in features}
    orientations = []
    for rank, oi in enumerate(order):
        fr = design_frames[oi]
        v = folded[fr]
        theta = float(np.degrees(np.arccos(np.clip(v[2], -1, 1))))
        phi = float(np.degrees(np.arctan2(v[1], v[0])))
        members = np.where(assign == oi)[0]
        orientations.append({
            "id": make_id(rank, v),
            "abc": [float(x) for x in v],           # body-frame field unit vector
            "theta_deg": round(theta, 3),           # polar from body +z
            "phi_deg": round(phi, 3),               # azimuth atan2(b,a)
            "feature": feat_of.get(fr),
            "representative_frame": int(fr),
            "dwell_fraction": float(weights[oi]),   # Stage-2 Voronoi occupancy
            "n_member_frames": int(members.size),
            "member_frames_sample": members[:: max(1, members.size // 8)][:8].tolist(),
        })

    manifest = {
        "meta": {
            "source_trajectory": args.trajectory,
            "n_frames": int(n),
            "timestep_au": TIMESTEP_AU,
            "timestep_fs": round(TIMESTEP_AU * AU_TIME_FS, 4),
            "lambda_magnitude": args.lambda_mag,
            "lab_field_hat": [float(x) for x in f_hat],
            "annotation_residual_deg": round(float(resid), 5),
            "intermediates": args.intermediates,
            "n_orientations": len(orientations),
            "n_cells": len(orientations) * len(args.intermediates),
            "pc1_variance_fraction": round(float(var[0]), 3),
            "frame_convention": (
                "orthonormal body frame: x=ipso_C->N; z=ring-normal (first-3-C "
                "triplet cross product) orthogonalized to x, molecular sign fixed "
                "by carbon ordering (NO lab flip); y=z cross x."),
            "symmetry_fold": (
                "site mirror sigma_v folded (b->|b|): one intermediate at "
                "(a,|b|,c) represents both mirror orientations. Global parity "
                "E(lambda)=E(-lambda) unused (opposite face never visited)."),
            "placement": (
                "STAGE 1 dwell-BLIND: farthest-point coverage of the support, "
                "with one design point pinned at each PES feature (2 basins, 2 "
                "major hills, 1 minor hill). Dwell re-enters only at STAGE 2 "
                "reweighting (see frames.npz)."),
            "campaign_contract": (
                "A cell is (intermediate, orientation.id). Build lambda per cell "
                "with orientation_transfer.lambda_for_intermediate(orientation "
                "['abc'], symbols, intermediate_ref_coords, lambda_magnitude). "
                "That helper asserts Br-on-minus-z and frame orthonormality; a "
                "raised AssertionError means the intermediate atom ordering "
                "disagrees with the reactant frame -- fix before spending CPU."),
        },
        "orientations": orientations,
    }
    with open(f"{args.out_dir}/manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[write] manifest.json ({len(orientations)} orientations, "
          f"{manifest['meta']['n_cells']} cells)")

    # convenience: cell id list for campaign --list-ids style consumption
    with open(f"{args.out_dir}/cell_ids.txt", "w") as fh:
        for o in orientations:
            for inter in args.intermediates:
                fh.write(f"{inter}__{o['id']}\n")

    _verification_figure(args.out_dir, folded, s, E_kcal, cent, Es,
                         design_frames, feat_frames, assign)
    print(f"[write] manifest_check.png, cell_ids.txt, frames.npz")


def _verification_figure(out, folded, s, E_kcal, cent, Es, design, feats, assign):
    def lambert(p):
        px, py, pz = p[:, 0], p[:, 1], p[:, 2]
        k = np.sqrt(2.0 / (1.0 - pz + 1e-12))       # about -z pole
        return k * px, k * py

    fig, ax = plt.subplots(1, 2, figsize=(15, 6.2))

    X, Y = lambert(folded)
    ax[0].hexbin(X, Y, gridsize=45, bins="log", cmap="Blues",
                 extent=(-1.5, 1.5, -0.2, 1.6))
    Xd, Yd = X[design], Y[design]
    ax[0].scatter(Xd, Yd, s=70, facecolor="none", edgecolor="k", lw=1.4,
                  label="design orientation")
    ax[0].scatter(X[feats], Y[feats], s=130, marker="*", color="crimson",
                  zorder=5, label="pinned PES feature")
    th = np.linspace(0, np.pi, 100)
    ax[0].plot(np.sqrt(2) * np.cos(th), np.sqrt(2) * np.sin(th), "k-", lw=0.6)
    ax[0].set_title("Design orientations over dwell density\n"
                    "(folded to b>=0, Lambert equal-area)")
    ax[0].set_xlabel("body x"); ax[0].set_ylabel("body |y|")
    ax[0].set_aspect("equal"); ax[0].legend(loc="upper right", fontsize=9)

    ax[1].plot(cent, Es, "k-", lw=2, label="orientational PES E(s)")
    sd = s[design]; Ed = np.interp(sd, cent, Es)
    ax[1].scatter(sd, Ed, s=60, facecolor="none", edgecolor="steelblue", lw=1.6,
                  zorder=4, label="design points")
    sf = s[feats]; Ef = np.interp(sf, cent, Es)
    ax[1].scatter(sf, Ef, s=150, marker="*", color="crimson", zorder=5,
                  label="pinned features")
    ax[1].set_xlabel("path coordinate s (PC1)")
    ax[1].set_ylabel("E - E_min (kcal/mol)")
    ax[1].set_title("Coverage along the swing coordinate\n"
                    "features pinned, gaps filled by farthest-point")
    ax[1].legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{out}/manifest_check.png", dpi=130)


if __name__ == "__main__":
    main()
