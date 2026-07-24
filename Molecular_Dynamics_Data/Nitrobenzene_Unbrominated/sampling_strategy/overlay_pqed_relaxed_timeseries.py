#!/usr/bin/env python3
"""Sandbox overlay of pQED surfaces with relaxed sampled-orientation energies.

This intentionally includes stalled relaxed rows when they have
final_energy_hartree, marking them as provisional. ZPE columns are present in
the source summary but are not included here.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/private/tmp") / "mplconfig-codex"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator


HA_TO_KCAL = 627.5094740631
TS_TO_FS = 25.0 * 0.02418884

DATA_ROOT = Path("/Users/jfoley19/Code/nbr_data/Molecular_Dynamics_Data")
PQED_DIR = DATA_ROOT / "Nitrobenzene_Unbrominated" / "QED_DFT_wb97x_d"
TIMESERIES_SCRIPT = PQED_DIR / "plot_timeseries_with_dwell_times.py"
PQED_GRID = PQED_DIR / "isomer_Nel_49_Nph_10_total_energies.dat"
RELAXED_CSV = DATA_ROOT / "Direction_A_MD_Sampled_Orientations" / "qed_dft_energy_summary.csv"
OUT_DIR = Path(__file__).resolve().parent


def fold_theta(theta: np.ndarray | pd.Series) -> np.ndarray:
    theta = np.asarray(theta, float)
    return np.where(theta > 100.0, 180.0 - theta, theta)


def parse_md(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("Step"):
            parts = line.split()
            theta_raw = float(parts[4].split("=")[1])
            rows.append(
                {
                    "step": int(parts[1]),
                    "time_fs": int(parts[1]) * TS_TO_FS,
                    "e_md": float(parts[2].split("=")[1]),
                    "phi": float(parts[3].split("=")[1]),
                    "theta_raw": theta_raw,
                    "theta": float(fold_theta([theta_raw])[0]),
                }
            )
    return pd.DataFrame(rows)


def current_md_file_from_timeseries_script() -> Path:
    text = TIMESERIES_SCRIPT.read_text()
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("MD_FILE"):
            _, rhs = stripped.split("=", 1)
            return PQED_DIR / rhs.split("#", 1)[0].strip().strip("\"'")
    raise RuntimeError(f"Could not find MD_FILE in {TIMESERIES_SCRIPT}")


def load_pqed_grid(path: Path) -> pd.DataFrame:
    grid = pd.read_csv(path, sep=r"\s+", skiprows=[1])
    for col in ["theta", "phi", "Para_E", "Ortho_E", "Meta_E"]:
        grid[col] = pd.to_numeric(grid[col], errors="coerce")
    grid["pqed_dE_om"] = (grid["Ortho_E"] - grid["Meta_E"]) * HA_TO_KCAL
    grid["pqed_dE_pm"] = (grid["Para_E"] - grid["Meta_E"]) * HA_TO_KCAL
    return grid


def add_pqed_to_md(grid: pd.DataFrame, md: pd.DataFrame) -> pd.DataFrame:
    theta_vals = np.sort(grid["theta"].unique())
    phi_vals = np.sort(grid["phi"].unique())

    def interp(col: str) -> RegularGridInterpolator:
        values = grid.pivot(index="theta", columns="phi", values=col).values
        return RegularGridInterpolator((theta_vals, phi_vals), values, bounds_error=False, fill_value=None)

    pts = md[["theta", "phi"]].to_numpy(float)
    out = md.copy()
    out["pqed_dE_om"] = interp("pqed_dE_om")(pts)
    out["pqed_dE_pm"] = interp("pqed_dE_pm")(pts)
    return out


def load_relaxed_pairs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df["final_energy_hartree"].notna()].copy()
    df["theta_plot"] = fold_theta(df["theta"])
    df["pair_quality"] = np.where(
        df["has_optimized_xyz"].eq(True) & df["converged"].eq(True),
        "converged",
        "provisional_stalled",
    )
    pivot_e = df.pivot_table(index=["theta", "phi", "theta_plot"], columns="isomer", values="final_energy_hartree", aggfunc="first")
    pivot_q = df.pivot_table(index=["theta", "phi", "theta_plot"], columns="isomer", values="pair_quality", aggfunc="first")
    out = pivot_e.reset_index()
    out.columns.name = None
    q = pivot_q.reset_index()
    q.columns.name = None
    for isomer in ["ortho", "meta", "para"]:
        if isomer not in out:
            out[isomer] = np.nan
        if isomer not in q:
            q[isomer] = np.nan
    out["relaxed_dE_om"] = (out["ortho"] - out["meta"]) * HA_TO_KCAL
    out["relaxed_dE_pm"] = (out["para"] - out["meta"]) * HA_TO_KCAL
    out["quality_om"] = np.where((q["ortho"] == "converged") & (q["meta"] == "converged"), "converged", "provisional_stalled")
    out["quality_pm"] = np.where((q["para"] == "converged") & (q["meta"] == "converged"), "converged", "provisional_stalled")
    return out


def fourier_features(theta_deg: np.ndarray, phi_deg: np.ndarray, order: int = 2) -> np.ndarray:
    theta = np.radians(theta_deg)
    phi = np.radians(phi_deg)
    cols = [np.ones_like(theta)]
    for n in range(1, order + 1):
        cols.extend([np.cos(n * theta), np.sin(n * theta)])
    for m in range(1, order + 1):
        cols.extend([np.cos(m * phi), np.sin(m * phi)])
    for n in range(1, order + 1):
        for m in range(1, order + 1):
            cols.extend(
                [
                    np.cos(n * theta) * np.cos(m * phi),
                    np.cos(n * theta) * np.sin(m * phi),
                    np.sin(n * theta) * np.cos(m * phi),
                    np.sin(n * theta) * np.sin(m * phi),
                ]
            )
    return np.column_stack(cols)


def fit_fourier(points: pd.DataFrame, value_col: str, alpha: float = 1e-2):
    fit_df = points[points[value_col].notna()].copy()
    x = fourier_features(fit_df["theta_plot"].to_numpy(float), fit_df["phi"].to_numpy(float))
    y = fit_df[value_col].to_numpy(float)
    scale = x.std(axis=0)
    scale[0] = 1.0
    scale[scale == 0] = 1.0
    xs = x / scale
    penalty = np.eye(xs.shape[1])
    penalty[0, 0] = 0.0
    beta = np.linalg.solve(xs.T @ xs + alpha * penalty, xs.T @ y)

    def predict(theta_deg, phi_deg):
        xx = fourier_features(np.asarray(theta_deg, float), np.asarray(phi_deg, float)) / scale
        return xx @ beta

    pred_train = predict(fit_df["theta_plot"], fit_df["phi"])
    rmse = float(np.sqrt(np.mean((pred_train - y) ** 2))) if len(y) else np.nan
    return predict, {"n_points": int(len(y)), "train_rmse_kcal_mol": rmse, "alpha": alpha}


def nearest_times(md: pd.DataFrame, points: pd.DataFrame) -> pd.DataFrame:
    out = points.copy()
    md_theta = md["theta"].to_numpy(float)
    md_phi = md["phi"].to_numpy(float)
    nearest_step, nearest_time, nearest_dist = [], [], []
    for _, row in out.iterrows():
        dphi = np.abs(md_phi - row["phi"])
        dphi = np.minimum(dphi, 360.0 - dphi)
        dist = np.hypot(md_theta - row["theta_plot"], dphi)
        idx = int(np.argmin(dist))
        nearest_step.append(int(md.iloc[idx]["step"]))
        nearest_time.append(float(md.iloc[idx]["time_fs"]))
        nearest_dist.append(float(dist[idx]))
    out["nearest_trajectory_step"] = nearest_step
    out["nearest_trajectory_time_fs"] = nearest_time
    out["nearest_trajectory_angular_dist_deg"] = nearest_dist
    return out


def surface_overlay(grid: pd.DataFrame, points: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), constrained_layout=True)
    panels = [
        ("pQED ortho - meta with relaxed samples", "pqed_dE_om", "relaxed_dE_om", "quality_om"),
        ("pQED para - meta with relaxed samples", "pqed_dE_pm", "relaxed_dE_pm", "quality_pm"),
    ]
    vlim = math.ceil(float(np.nanmax(np.abs(grid[["pqed_dE_om", "pqed_dE_pm"]].to_numpy()))) * 10) / 10
    for ax, (title, gcol, pcol, qcol) in zip(axes, panels):
        piv = grid.pivot(index="theta", columns="phi", values=gcol)
        im = ax.pcolormesh(piv.columns, piv.index, piv.values, cmap="coolwarm", vmin=-vlim, vmax=vlim, shading="auto")
        ax.contour(piv.columns, piv.index, piv.values, levels=[-5, 0, 5], colors="black", linewidths=[0.7, 1.0, 0.7], alpha=0.55)
        have = points[pcol].notna()
        conv = have & points[qcol].eq("converged")
        prov = have & ~points[qcol].eq("converged")
        ax.scatter(points.loc[prov, "phi"], points.loc[prov, "theta_plot"], s=82, facecolors="none", edgecolors="black", linewidths=1.3, label="provisional")
        ax.scatter(points.loc[conv, "phi"], points.loc[conv, "theta_plot"], s=72, c=points.loc[conv, pcol], cmap="coolwarm", vmin=-vlim, vmax=vlim, edgecolors="black", linewidths=0.9, label="converged")
        ax.scatter(points.loc[prov, "phi"], points.loc[prov, "theta_plot"], s=35, c=points.loc[prov, pcol], cmap="coolwarm", vmin=-vlim, vmax=vlim, edgecolors="none")
        ax.set_title(title)
        ax.set_xlabel("phi, deg")
        ax.set_ylabel("folded theta, deg")
        ax.set_xlim(0, 130)
        ax.set_ylim(0, 85)
        ax.legend(loc="upper right", fontsize=9)
    cbar = fig.colorbar(im, ax=axes, shrink=0.86)
    cbar.set_label("Delta E, kcal/mol")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def fitted_surface(points: pd.DataFrame, fit_om, fit_pm, out_path: Path) -> None:
    phi = np.linspace(0, 130, 150)
    theta = np.linspace(0, 85, 120)
    pp, tt = np.meshgrid(phi, theta)
    z_om = fit_om(tt.ravel(), pp.ravel()).reshape(tt.shape)
    z_pm = fit_pm(tt.ravel(), pp.ravel()).reshape(tt.shape)
    vlim = math.ceil(float(np.nanmax(np.abs([z_om, z_pm]))) * 10) / 10
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), constrained_layout=True)
    panels = [
        ("Fourier fit: relaxed ortho - meta", z_om, "relaxed_dE_om", "quality_om"),
        ("Fourier fit: relaxed para - meta", z_pm, "relaxed_dE_pm", "quality_pm"),
    ]
    for ax, (title, z, pcol, qcol) in zip(axes, panels):
        im = ax.pcolormesh(phi, theta, z, cmap="coolwarm", vmin=-vlim, vmax=vlim, shading="auto")
        ax.contour(phi, theta, z, levels=[-5, 0, 5], colors="black", linewidths=[0.7, 1.0, 0.7], alpha=0.55)
        have = points[pcol].notna()
        conv = have & points[qcol].eq("converged")
        prov = have & ~points[qcol].eq("converged")
        ax.scatter(points.loc[prov, "phi"], points.loc[prov, "theta_plot"], s=78, facecolors="none", edgecolors="black", linewidths=1.25)
        ax.scatter(points.loc[conv, "phi"], points.loc[conv, "theta_plot"], s=72, c=points.loc[conv, pcol], cmap="coolwarm", vmin=-vlim, vmax=vlim, edgecolors="black", linewidths=0.9)
        ax.scatter(points.loc[prov, "phi"], points.loc[prov, "theta_plot"], s=32, c=points.loc[prov, pcol], cmap="coolwarm", vmin=-vlim, vmax=vlim, edgecolors="none")
        ax.set_title(title)
        ax.set_xlabel("phi, deg")
        ax.set_ylabel("folded theta, deg")
    cbar = fig.colorbar(im, ax=axes, shrink=0.86)
    cbar.set_label("Delta E, kcal/mol")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def timeseries_overlay(md: pd.DataFrame, points: pd.DataFrame, fit_om, fit_pm, out_path: Path) -> None:
    ts = md.iloc[10:].copy()
    ts["fit_dE_om"] = fit_om(ts["theta"], ts["phi"])
    ts["fit_dE_pm"] = fit_pm(ts["theta"], ts["phi"])
    near = nearest_times(md, points)
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    ax.plot(ts["time_fs"], ts["pqed_dE_om"], color="#0072B2", lw=2.0, label="pQED ortho - meta")
    ax.plot(ts["time_fs"], ts["pqed_dE_pm"], color="#E69F00", lw=2.0, ls="--", label="pQED para - meta")
    ax.plot(ts["time_fs"], ts["fit_dE_om"], color="#0072B2", lw=1.3, alpha=0.55, ls=":", label="relaxed Fourier fit ortho - meta")
    ax.plot(ts["time_fs"], ts["fit_dE_pm"], color="#E69F00", lw=1.3, alpha=0.55, ls=":", label="relaxed Fourier fit para - meta")
    for ycol, qcol, color, marker in [
        ("relaxed_dE_om", "quality_om", "#0072B2", "o"),
        ("relaxed_dE_pm", "quality_pm", "#E69F00", "s"),
    ]:
        have = near[ycol].notna()
        conv = have & near[qcol].eq("converged")
        prov = have & ~near[qcol].eq("converged")
        ax.scatter(near.loc[prov, "nearest_trajectory_time_fs"], near.loc[prov, ycol], marker=marker, s=62, facecolors="none", edgecolors=color, linewidths=1.3)
        ax.scatter(near.loc[conv, "nearest_trajectory_time_fs"], near.loc[conv, ycol], marker=marker, s=54, c=color, edgecolors="black", linewidths=0.6)
    ax.axhline(0.0, color="black", lw=1.0, alpha=0.7)
    ax.axhline(5.0, color="gray", ls=":", lw=0.9)
    ax.axhline(-5.0, color="gray", ls=":", lw=0.9)
    ax.set_xlabel("Trajectory time, fs")
    ax.set_ylabel("Delta E, kcal/mol")
    ax.set_title("pQED timeseries with provisional relaxed-data overlay")
    ax.grid(True, ls=":", alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return near


def main() -> None:
    md_file = current_md_file_from_timeseries_script()
    grid = load_pqed_grid(PQED_GRID)
    md = add_pqed_to_md(grid, parse_md(md_file))
    points = load_relaxed_pairs(RELAXED_CSV)
    fit_om, stats_om = fit_fourier(points, "relaxed_dE_om")
    fit_pm, stats_pm = fit_fourier(points, "relaxed_dE_pm")

    surface_overlay(grid, points, OUT_DIR / "pqed_surface_with_relaxed_samples.png")
    fitted_surface(points, fit_om, fit_pm, OUT_DIR / "relaxed_fourier_fit_surfaces.png")
    near = timeseries_overlay(md, points, fit_om, fit_pm, OUT_DIR / "trajectory_pqed_relaxed_timeseries_overlay.png")

    near.to_csv(OUT_DIR / "relaxed_points_nearest_trajectory_times.csv", index=False)
    points.to_csv(OUT_DIR / "relaxed_sandbox_pair_points.csv", index=False)
    summary = {
        "notes": [
            "Sandbox includes stalled rows with final_energy_hartree as provisional relaxed electronic energies.",
            "Converged/optimized pairs are still marked separately.",
            "ZPE columns are not included in this pass.",
            "Theta values are folded with the same theta > 100 -> 180 - theta convention used by plot_timeseries_with_dwell_times.py.",
            "Fourier fit is a low-order ridge-regularized sinusoidal fit and should be treated as provisional while meta remains sparse.",
        ],
        "md_file": str(md_file),
        "counts": {
            "relaxed_om_pairs_all_available": int(points["relaxed_dE_om"].notna().sum()),
            "relaxed_pm_pairs_all_available": int(points["relaxed_dE_pm"].notna().sum()),
            "relaxed_om_pairs_converged": int((points["relaxed_dE_om"].notna() & points["quality_om"].eq("converged")).sum()),
            "relaxed_pm_pairs_converged": int((points["relaxed_dE_pm"].notna() & points["quality_pm"].eq("converged")).sum()),
        },
        "fit_stats": {
            "ortho_meta": stats_om,
            "para_meta": stats_pm,
        },
        "nearest_trajectory_distance_deg": {
            "min": float(near["nearest_trajectory_angular_dist_deg"].min()),
            "median": float(near["nearest_trajectory_angular_dist_deg"].median()),
            "max": float(near["nearest_trajectory_angular_dist_deg"].max()),
        },
    }
    (OUT_DIR / "pqed_relaxed_overlay_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
