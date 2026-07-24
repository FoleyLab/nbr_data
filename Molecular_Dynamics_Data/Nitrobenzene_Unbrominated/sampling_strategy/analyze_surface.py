#!/usr/bin/env python3
"""
analyze_surface.py -- Stage-2 analysis for the MD orientation sweep.
====================================================================

Reads the completed cells in runs/, fits a smooth orientation surface E+ZPE for
each intermediate, and answers the scientific question:

    dE_om(orientation) = (E+ZPE)_ortho - (E+ZPE)_meta          (and vs para)

    "For each field orientation the reactant visits, which Wheland intermediate
     is most energetically accessible?"

It also (a) reweights that surface by the trajectory's true dwell occupancy from
frames.npz -- the Stage-2 step that was deliberately kept OUT of sampling -- and
(b) emits a batch-2 acquisition score to place the next ~N orientations where the
answer is least certain and matters most.

Design points were placed dwell-BLIND for coverage; dwell enters ONLY here. The
same fitted surface can be reweighted under any ensemble (microcanonical dwell,
Boltzmann, or flat) without recomputing anything.

Method
------
Per intermediate, fit a Gaussian process on the 3D unit vectors (a,b,c) with an
RBF kernel + white noise. On a small spherical patch the chordal RBF is an
excellent proxy for a great-circle kernel and keeps us on scikit-learn's tested
path. Predictions and 1-sigma bands transfer to dE by summing variances (the
GP fits are independent per intermediate).

Acquisition for batch 2 (per candidate frame orientation x):
    score(x) = sigma_dE(x) * dwell_density(x) * boundary(x)
where boundary(x) = exp(-(mu_dE(x)/tau)^2) up-weights the dE=0 flip contour --
the only place a sign error changes the scientific conclusion.

Usage
-----
    python analyze_surface.py [--runs runs] [--frames frames.npz]
        [--manifest manifest.json] [--n-batch2 15] [--out-dir .]
        [--boltzmann-T 0]        # 0 => microcanonical dwell weighting
"""

import argparse, json, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HA2KCAL = 627.5094740631


# --------------------------------------------------------------------------- loading
def load_cells(runs_dir):
    """Collect completed cells: {intermediate: {ori_id: (abc, E+ZPE, n_imag)}}."""
    out = {}
    for done in glob.glob(os.path.join(runs_dir, "*", "DONE")):
        cell_dir = os.path.dirname(done)
        cell_id = os.path.basename(cell_dir)
        inter, ori_id = cell_id.split("__", 1)
        try:
            freq = json.load(open(os.path.join(cell_dir, "frequencies.json")))
            opt = json.load(open(os.path.join(cell_dir, "opt_status.json")))
        except FileNotFoundError:
            continue
        out.setdefault(inter, {})[ori_id] = {
            "abc": np.array(opt["abc"]),
            "e_zpe": freq.get("e_plus_zpe_hartree",
                              freq["final_energy_hartree"] + freq["zpe_hartree"]),
            "n_imag": freq.get("n_imaginary", 0),
        }
    return out


# ----------------------------------------------------------------------------- GP fit
def fit_gp(X, y):
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
    y0 = y.mean()
    kernel = (ConstantKernel(1.0, (1e-6, 1e6)) * RBF(0.5, (0.02, 50.0))
              + WhiteKernel(1e-6, (1e-12, 1e-2)))
    gp = GaussianProcessRegressor(kernel=kernel, normalize_y=False,
                                  n_restarts_optimizer=4, alpha=1e-10)
    gp.fit(X, y - y0)
    return gp, y0


def predict(gp, y0, X):
    mu, sd = gp.predict(X, return_std=True)
    return mu + y0, sd


# -------------------------------------------------------------------- reweight helper
def ensemble_weights(E_hartree, dwell, T):
    """T<=0 => microcanonical dwell; T>0 => Boltzmann in dwell-sampled ensemble."""
    if T is None or T <= 0:
        return dwell / dwell.sum()
    kT = 3.166808e-6 * T                     # Hartree per Kelvin
    w = dwell * np.exp(-(E_hartree - E_hartree.min()) / kT)
    return w / w.sum()


# ---------------------------------------------------------------------------- reports
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--frames", default="frames.npz")
    ap.add_argument("--manifest", default="manifest.json")
    ap.add_argument("--n-batch2", type=int, default=15)
    ap.add_argument("--boltzmann-T", type=float, default=0.0)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    cells = load_cells(args.runs)
    if not cells:
        print(f"no completed cells in {args.runs}/ -- run campaign_md.py first.")
        return
    inters = sorted(cells)
    print("intermediates with data:", {k: len(v) for k, v in cells.items()})

    frames = np.load(args.frames)
    fabc, dwell, E_traj = frames["abc"], frames["dwell"], frames["E_hartree"]
    fold = fabc.copy(); fold[:, 1] = np.abs(fold[:, 1])       # match manifest fold

    # ---- fit a GP per intermediate on its computed orientations
    gps = {}
    for inter in inters:
        d = cells[inter]
        X = np.array([v["abc"] for v in d.values()])
        y = np.array([v["e_zpe"] for v in d.values()])
        n_imag_tot = sum(v["n_imag"] for v in d.values())
        gp, y0 = fit_gp(X, y)
        gps[inter] = (gp, y0)
        print(f"  {inter}: {len(y)} cells, {n_imag_tot} with imaginary modes, "
              f"kernel={gp.kernel_}")

    # ---- CSV: predicted E+ZPE and dE at every trajectory frame (Hartree)
    preds = {inter: predict(gp, y0, fold) for inter, (gp, y0) in gps.items()}
    csv = os.path.join(args.out_dir, "surface_per_frame.csv")
    with open(csv, "w") as fh:
        cols = sum(([f"E_{i}", f"sd_{i}"] for i in inters), [])
        fh.write("frame,a,b,c,dwell," + ",".join(cols))
        if {"ortho", "meta"} <= set(inters):
            fh.write(",dE_ortho_meta,sd_dE")
        fh.write("\n")
        mu_o, sd_o = preds.get("ortho", (None, None))
        mu_m, sd_m = preds.get("meta", (None, None))
        for t in range(len(fabc)):
            row = [t, *fabc[t], dwell[t]]
            for i in inters:
                row += [preds[i][0][t], preds[i][1][t]]
            line = ",".join(f"{x:.6g}" for x in row)
            if mu_o is not None and mu_m is not None:
                dE = (mu_o[t] - mu_m[t]) * HA2KCAL
                sdE = np.hypot(sd_o[t], sd_m[t]) * HA2KCAL
                line += f",{dE:.6g},{sdE:.6g}"
            fh.write(line + "\n")
    print(f"[write] {csv}")

    # ---- dwell-weighted preference summary
    if {"ortho", "meta"} <= set(inters):
        w = ensemble_weights(E_traj, dwell, args.boltzmann_T)
        dE = (mu_o - mu_m) * HA2KCAL
        frac_ortho = float(w[dE < 0].sum())
        ens = "Boltzmann T=%g K" % args.boltzmann_T if args.boltzmann_T > 0 else "microcanonical dwell"
        print(f"\n[{ens}] dwell-weighted <dE_(ortho-meta)> = "
              f"{np.average(dE, weights=w):+.3f} kcal/mol; "
              f"ortho favored in {100*frac_ortho:.1f}% of occupancy")

        # ---- batch-2 acquisition
        tau = 0.5                                  # kcal/mol; width of flip band
        # dwell density at each frame ~ its own weight (already per-frame)
        boundary = np.exp(-(dE / tau) ** 2)
        score = np.hypot(sd_o, sd_m) * HA2KCAL * (w / w.max()) * (0.2 + boundary)
        # greedily pick spread-out high-score frames (min great-circle sep)
        picks, chosen = [], np.ones(len(fabc), bool)
        idx_sorted = np.argsort(-score)
        for t in idx_sorted:
            if len(picks) >= args.n_batch2:
                break
            if all(np.arccos(np.clip(fold[t] @ fold[p], -1, 1)) > 0.12 for p in picks):
                picks.append(int(t))
        json.dump([{"frame": p, "abc": fold[p].tolist(),
                    "score": float(score[p])} for p in picks],
                  open(os.path.join(args.out_dir, "batch2_candidates.json"), "w"),
                  indent=2)
        print(f"[write] batch2_candidates.json ({len(picks)} orientations)")

        _figure(args.out_dir, fold, dE, np.hypot(sd_o, sd_m) * HA2KCAL, w,
                np.array([fold[p] for p in picks]))


def _figure(out, fold, dE, sdE, w, batch2):
    def lambert(p):
        k = np.sqrt(2.0 / (1.0 - p[:, 2] + 1e-12))
        return k * p[:, 0], k * p[:, 1]
    X, Y = lambert(fold)
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    sc = ax[0].scatter(X, Y, c=dE, s=8, cmap="coolwarm",
                       vmin=-np.abs(dE).max(), vmax=np.abs(dE).max())
    plt.colorbar(sc, ax=ax[0], label="dE (ortho-meta), kcal/mol")
    # mark predicted flip contour: frames where |dE| small
    flip = np.abs(dE) < 0.25
    ax[0].scatter(X[flip], Y[flip], s=14, facecolor="none", edgecolor="k", lw=0.4)
    if len(batch2):
        Xb, Yb = lambert(batch2)
        ax[0].scatter(Xb, Yb, marker="^", s=90, color="lime", edgecolor="k",
                      zorder=5, label="batch-2 pick")
        ax[0].legend()
    ax[0].set_title("dE(ortho-meta) over visited orientations\n"
                    "black ring = predicted flip contour")
    ax[0].set_aspect("equal"); ax[0].set_xlabel("body x"); ax[0].set_ylabel("body |y|")

    order = np.argsort(w)[::-1]
    ax[1].scatter(dE, sdE, c=w, s=10, cmap="viridis")
    ax[1].axvline(0, color="k", lw=0.7)
    ax[1].set_xlabel("dE (ortho-meta), kcal/mol")
    ax[1].set_ylabel("GP 1-sigma on dE, kcal/mol")
    ax[1].set_title("Uncertainty vs preference\n(color = dwell weight)")
    plt.tight_layout()
    plt.savefig(os.path.join(out, "surface_analysis.png"), dpi=130)
    print(f"[write] surface_analysis.png")


if __name__ == "__main__":
    main()
