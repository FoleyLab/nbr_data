#!/usr/bin/env python3
"""Plot relaxed and unrelaxed sampled-orientation energy differences."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HA_TO_KCAL = 627.5094740631

DATA_ROOT = Path("/Users/jfoley19/Code/nbr_data/Molecular_Dynamics_Data")
SAMPLING = DATA_ROOT / "Nitrobenzene_Unbrominated" / "sampling_strategy"
RELAXED_CSV = (
    DATA_ROOT
    / "Direction_A_MD_Sampled_Orientations"
    / "qed_dft_energy_summary.csv"
)
UNRELAXED_CSV = (
    DATA_ROOT
    / "Direction_A_MD_Sampled_Orientations"
    / "UNRELAXED_Sampled_Orientations"
    / "direction_A_unrelaxed_results"
    / "direction_A_unrelaxed_energies.csv"
)


def lambert(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Same simple folded-support projection used by the campaign scripts."""
    denom = 1.0 - points[:, 2] + 1e-12
    k = np.sqrt(2.0 / denom)
    return k * points[:, 0], k * points[:, 1]


def load_manifest() -> tuple[dict, pd.DataFrame, dict[str, int]]:
    manifest = json.loads((SAMPLING / "manifest.json").read_text())
    rows = []
    for idx, ori in enumerate(manifest["orientations"]):
        rows.append(
            {
                "ori_index": idx,
                "ori_id": ori["id"],
                "a": ori["abc"][0],
                "b": ori["abc"][1],
                "c": ori["abc"][2],
                "theta_manifest": ori["theta_deg"],
                "phi_manifest": ori["phi_deg"],
                "representative_frame": ori["representative_frame"],
                "dwell_fraction": ori["dwell_fraction"],
                "n_member_frames": ori["n_member_frames"],
                "feature": ori.get("feature") or "",
            }
        )
    ori_df = pd.DataFrame(rows)
    return manifest, ori_df, {r["ori_id"]: int(r["ori_index"]) for r in rows}


def match_orientations(df: pd.DataFrame, ori_df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Attach manifest orientation IDs by matching actual lambda geometry.

    Earlier direction labels use a rounded angular convention where the x sign can
    differ from the folded manifest representation, so the robust match is by
    vector overlap against both the recorded vector and its x-reflected partner.
    """
    out = df.copy()
    manifest_vectors = ori_df[["a", "b", "c"]].to_numpy(float)
    ori_ids = ori_df["ori_id"].to_numpy()
    matches = []
    for _, row in out.iterrows():
        v = row[cols].to_numpy(float)
        v /= np.linalg.norm(v)
        candidates = [v, np.array([-v[0], v[1], v[2]])]
        best_dot, best_idx = -np.inf, -1
        for candidate in candidates:
            dots = manifest_vectors @ candidate
            idx = int(np.argmax(dots))
            if dots[idx] > best_dot:
                best_dot, best_idx = float(dots[idx]), idx
        matches.append((ori_ids[best_idx], best_idx, best_dot))
    out["ori_id"] = [m[0] for m in matches]
    out["ori_index"] = [m[1] for m in matches]
    out["orientation_match_dot"] = [m[2] for m in matches]
    return out


def pivot_energy(df: pd.DataFrame, energy_col: str, prefix: str) -> pd.DataFrame:
    piv = df.pivot_table(index="ori_id", columns="isomer", values=energy_col, aggfunc="first")
    piv = piv.rename(columns={c: f"{prefix}_E_{c}" for c in piv.columns})
    piv = piv.reset_index()
    for pair, left, right in [
        ("ortho_meta", f"{prefix}_E_ortho", f"{prefix}_E_meta"),
        ("para_meta", f"{prefix}_E_para", f"{prefix}_E_meta"),
    ]:
        if left in piv and right in piv:
            piv[f"{prefix}_dE_{pair}_kcal_mol"] = (piv[left] - piv[right]) * HA_TO_KCAL
        else:
            piv[f"{prefix}_dE_{pair}_kcal_mol"] = np.nan
    return piv


def base_axis(ax, frames_folded: np.ndarray, ori_df: pd.DataFrame) -> None:
    x, y = lambert(frames_folded)
    ax.scatter(x, y, s=4, c="0.82", alpha=0.45, linewidths=0, rasterized=True)
    ox, oy = lambert(ori_df[["a", "b", "c"]].to_numpy(float))
    ax.scatter(ox, oy, s=30, facecolors="none", edgecolors="0.25", linewidths=0.55)
    feat = ori_df["feature"].astype(bool).to_numpy()
    if feat.any():
        ax.scatter(ox[feat], oy[feat], s=95, facecolors="none", edgecolors="black", linewidths=1.2)
    ax.set_aspect("equal")
    ax.set_xlabel("body x")
    ax.set_ylabel("body |y|")


def plot_map(
    merged: pd.DataFrame,
    frames_folded: np.ndarray,
    cols: list[str],
    titles: list[str],
    output: Path,
    vlim: float,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8), constrained_layout=True)
    points = merged[["a", "b", "c"]].to_numpy(float)
    px, py = lambert(points)
    cmap = "coolwarm"
    for ax, col, title in zip(axes, cols, titles):
        base_axis(ax, frames_folded, merged)
        available = merged[col].notna().to_numpy()
        sc = ax.scatter(
            px[available],
            py[available],
            c=merged.loc[available, col],
            s=82,
            cmap=cmap,
            vmin=-vlim,
            vmax=vlim,
            edgecolors="black",
            linewidths=0.45,
            zorder=3,
        )
        missing = ~available
        if missing.any():
            ax.scatter(
                px[missing],
                py[missing],
                s=58,
                facecolors="none",
                edgecolors="0.55",
                linewidths=0.9,
                zorder=2,
            )
        ax.axhline(0, color="0.88", lw=0.6, zorder=0)
        ax.axvline(0, color="0.88", lw=0.6, zorder=0)
        ax.set_title(f"{title}\n{available.sum()}/25 pairs available")
    cbar = fig.colorbar(sc, ax=axes, shrink=0.86)
    cbar.set_label("Delta E relative to meta, kcal/mol")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_relaxed_vs_unrelaxed(merged: pd.DataFrame, output: Path, vlim: float) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), constrained_layout=True)
    pairs = [
        ("ortho - meta", "unrelaxed_dE_ortho_meta_kcal_mol", "relaxed_dE_ortho_meta_kcal_mol"),
        ("para - meta", "unrelaxed_dE_para_meta_kcal_mol", "relaxed_dE_para_meta_kcal_mol"),
    ]
    for ax, (title, xcol, ycol) in zip(axes, pairs):
        ok = merged[xcol].notna() & merged[ycol].notna()
        colors = np.where(merged.loc[ok, "feature"].astype(bool), "black", "0.35")
        ax.scatter(merged.loc[ok, xcol], merged.loc[ok, ycol], s=70, c=colors, alpha=0.88)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.axvline(0, color="0.5", lw=0.8)
        lim = max(vlim, float(np.nanmax(np.abs(merged[[xcol, ycol]].to_numpy()))))
        ax.plot([-lim, lim], [-lim, lim], color="0.35", lw=1.0, ls="--")
        ax.set_xlim(-lim * 1.08, lim * 1.08)
        ax.set_ylim(-lim * 1.08, lim * 1.08)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{title}\n{ok.sum()} paired relaxed points")
        ax.set_xlabel("unrelaxed Delta E, kcal/mol")
        ax.set_ylabel("relaxed Delta E, kcal/mol")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    out_dir = Path(__file__).resolve().parent
    _, ori_df, _ = load_manifest()
    frames = np.load(SAMPLING / "frames.npz", allow_pickle=False)
    frames_folded = frames["folded"]

    unrelaxed = pd.read_csv(UNRELAXED_CSV)
    unrelaxed = unrelaxed[unrelaxed["status"].eq("ok")].copy()
    unrelaxed = match_orientations(unrelaxed, ori_df, ["lambda_x", "lambda_y", "lambda_z"])

    relaxed = pd.read_csv(RELAXED_CSV)
    relaxed = relaxed[
        relaxed["has_optimized_xyz"].eq(True)
        & relaxed["converged"].eq(True)
        & relaxed["final_energy_hartree"].notna()
    ].copy()
    relaxed = match_orientations(relaxed, ori_df, ["Ex", "Ey", "Ez"])

    unrelaxed_piv = pivot_energy(unrelaxed, "energy_hartree", "unrelaxed")
    relaxed_piv = pivot_energy(relaxed, "final_energy_hartree", "relaxed")
    merged = ori_df.merge(unrelaxed_piv, on="ori_id", how="left").merge(relaxed_piv, on="ori_id", how="left")

    merged["relaxed_electronic_energy_available"] = merged[
        ["relaxed_E_ortho", "relaxed_E_meta", "relaxed_E_para"]
    ].notna().any(axis=1)
    merged.to_csv(out_dir / "orientation_energy_differences_merged.csv", index=False)

    diff_cols = [
        "unrelaxed_dE_ortho_meta_kcal_mol",
        "unrelaxed_dE_para_meta_kcal_mol",
        "relaxed_dE_ortho_meta_kcal_mol",
        "relaxed_dE_para_meta_kcal_mol",
    ]
    vlim = math.ceil(float(np.nanmax(np.abs(merged[diff_cols].to_numpy()))) * 10) / 10

    plot_map(
        merged,
        frames_folded,
        ["unrelaxed_dE_ortho_meta_kcal_mol", "unrelaxed_dE_para_meta_kcal_mol"],
        ["Unrelaxed ortho - meta", "Unrelaxed para - meta"],
        out_dir / "unrelaxed_energy_differences_on_orientation_map.png",
        vlim,
    )
    plot_map(
        merged,
        frames_folded,
        ["relaxed_dE_ortho_meta_kcal_mol", "relaxed_dE_para_meta_kcal_mol"],
        ["Relaxed ortho - meta", "Relaxed para - meta"],
        out_dir / "relaxed_energy_differences_on_orientation_map.png",
        vlim,
    )
    plot_relaxed_vs_unrelaxed(
        merged,
        out_dir / "relaxed_vs_unrelaxed_energy_differences.png",
        vlim,
    )

    summary = {
        "notes": [
            "Relaxed plots use final_energy_hartree only.",
            "ZPE columns are present in the source relaxed CSV but intentionally not included in this pass.",
            "Only converged rows with optimized XYZ files are included for relaxed differences.",
        ],
        "counts": {
            "unrelaxed_rows_ok": int(len(unrelaxed)),
            "relaxed_rows_converged_optimized": int(len(relaxed)),
            "unrelaxed_ortho_meta_pairs": int(merged["unrelaxed_dE_ortho_meta_kcal_mol"].notna().sum()),
            "unrelaxed_para_meta_pairs": int(merged["unrelaxed_dE_para_meta_kcal_mol"].notna().sum()),
            "relaxed_ortho_meta_pairs": int(merged["relaxed_dE_ortho_meta_kcal_mol"].notna().sum()),
            "relaxed_para_meta_pairs": int(merged["relaxed_dE_para_meta_kcal_mol"].notna().sum()),
        },
        "min_orientation_match_dot": {
            "unrelaxed": float(unrelaxed["orientation_match_dot"].min()),
            "relaxed": float(relaxed["orientation_match_dot"].min()),
        },
        "color_scale_abs_kcal_mol": vlim,
    }
    (out_dir / "orientation_energy_differences_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
