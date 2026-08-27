#!/usr/bin/env python3
"""Build merged QED-DFT optimized energy tables and symmetry-expanded maps."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

try:
    from PIL import Image, ImageDraw, ImageFont
except ModuleNotFoundError:
    Image = ImageDraw = ImageFont = None


AU_TO_KCAL = 627.509
ISOMERS = ("ortho", "meta", "para")
INPUT_TEMPLATE = "grid_campaign_no_freq_{isomer}_status.csv"
OUTPUT_CSV = "grid_campaign_no_freq_opt_energies.csv"


def coordinate_key(theta: float, phi: float) -> tuple[float, float]:
    return (round(theta, 8), round(phi, 8))


def read_isomer_table(data_dir: Path, isomer: str) -> dict[tuple[float, float], float]:
    input_path = data_dir / INPUT_TEMPLATE.format(isomer=isomer)
    values: dict[tuple[float, float], float] = {}

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required_columns = {"theta", "phi", "E_opt_hartree"}
        missing = required_columns.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{input_path} is missing columns: {sorted(missing)}")

        for row in reader:
            theta = float(row["theta"])
            phi = float(row["phi"])
            energy = float(row["E_opt_hartree"])
            key = coordinate_key(theta, phi)
            if key in values:
                raise ValueError(f"Duplicate {isomer} row for theta={theta}, phi={phi}")
            values[key] = energy

    return values


def build_merged_rows(data_dir: Path) -> list[dict[str, float]]:
    tables = {isomer: read_isomer_table(data_dir, isomer) for isomer in ISOMERS}
    coordinate_sets = {isomer: set(table) for isomer, table in tables.items()}
    common_coordinates = set.intersection(*coordinate_sets.values())

    for isomer, coordinates in coordinate_sets.items():
        missing = sorted(common_coordinates.symmetric_difference(coordinates))
        if missing:
            sample = ", ".join(f"({t:g}, {p:g})" for t, p in missing[:5])
            raise ValueError(
                f"{isomer} coordinates do not match the common grid; first mismatch(es): {sample}"
            )

    rows = []
    for theta, phi in sorted(common_coordinates):
        rows.append(
            {
                "theta": theta,
                "phi": phi,
                "E_opt_ortho": tables["ortho"][(theta, phi)],
                "E_opt_meta": tables["meta"][(theta, phi)],
                "E_opt_para": tables["para"][(theta, phi)],
            }
        )
    return rows


def write_merged_csv(rows: list[dict[str, float]], output_path: Path) -> None:
    fieldnames = ["theta", "phi", "E_opt_ortho", "E_opt_meta", "E_opt_para"]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_lookup(rows: list[dict[str, float]]) -> dict[tuple[float, float], tuple[float, float, float]]:
    lookup = {}
    for row in rows:
        key = coordinate_key(row["theta"], row["phi"])
        lookup[key] = (row["E_opt_ortho"], row["E_opt_meta"], row["E_opt_para"])
    return lookup


def fetch_with_inversion_symmetry(
    theta: float,
    phi: float,
    lookup: dict[tuple[float, float], tuple[float, float, float]],
    theta0_data: tuple[float, float, float],
) -> tuple[float, float, float]:
    if np.isclose(theta, 0.0) or np.isclose(theta, 180.0):
        return theta0_data

    phi_wrapped = 0.0 if np.isclose(phi, 360.0) else phi
    direct_key = coordinate_key(theta, phi_wrapped)
    if direct_key in lookup:
        return lookup[direct_key]

    theta_inv = 180.0 - theta
    phi_inv = (phi_wrapped + 180.0) % 360.0
    inverted_key = coordinate_key(theta_inv, phi_inv)
    if inverted_key in lookup:
        return lookup[inverted_key]

    raise ValueError(
        f"No direct or inversion-symmetry value for theta={theta:g}, phi={phi:g}; "
        f"checked theta={theta_inv:g}, phi={phi_inv:g}"
    )


def symmetry_expanded_grids(
    rows: list[dict[str, float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    theta_values = sorted({row["theta"] for row in rows})
    phi_values = sorted({row["phi"] for row in rows})
    d_theta = min(np.diff(theta_values))
    d_phi = min(np.diff(phi_values))
    full_theta = np.arange(0.0, 180.0 + d_theta / 2.0, d_theta)
    full_phi = np.arange(0.0, 360.0 + d_phi / 2.0, d_phi)

    lookup = build_lookup(rows)
    theta0_data = lookup[coordinate_key(0.0, 0.0)]
    energy_grid = np.empty((len(full_theta), len(full_phi), 3), dtype=float)

    for i, theta in enumerate(full_theta):
        for j, phi in enumerate(full_phi):
            energy_grid[i, j, :] = fetch_with_inversion_symmetry(theta, phi, lookup, theta0_data)

    phi_grid, theta_grid = np.meshgrid(full_phi, full_theta)
    return theta_grid, phi_grid, energy_grid[:, :, 0], energy_grid[:, :, 1], energy_grid[:, :, 2]


def create_difference_map(
    theta_grid: np.ndarray,
    phi_grid: np.ndarray,
    difference_grid: np.ndarray,
    title: str,
    fallback_title: str,
    output_path: Path,
) -> None:
    if plt is None:
        create_difference_map_with_pillow(
            theta_grid,
            phi_grid,
            difference_grid,
            fallback_title,
            output_path,
        )
        return

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    vmax = float(np.nanmax(np.abs(difference_grid)))
    image = ax.pcolormesh(
        phi_grid,
        theta_grid,
        difference_grid,
        shading="gouraud",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )

    ax.set_title(title, fontsize=20, fontweight="bold", pad=20)
    ax.set_xlabel(r"Azimuthal Angle $\phi$ (deg.)", fontsize=18)
    ax.set_ylabel(r"Polar Angle $\theta$ (deg.)", fontsize=18)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 45, 90, 135, 180])
    ax.tick_params(axis="both", labelsize=16)
    ax.grid(True, linestyle="--", alpha=0.3)

    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Energy Difference (kcal/mol)", fontsize=18)
    colorbar.ax.tick_params(labelsize=16)

    fig.savefig(output_path, dpi=350)
    plt.close(fig)


def diverging_color(value: float, vmax: float) -> tuple[int, int, int]:
    if vmax <= 0.0:
        amount = 0.0
    else:
        amount = max(-1.0, min(1.0, value / vmax))

    white = np.array([247, 247, 247], dtype=float)
    blue = np.array([49, 130, 189], dtype=float)
    red = np.array([203, 24, 29], dtype=float)
    endpoint = red if amount >= 0.0 else blue
    color = white * (1.0 - abs(amount)) + endpoint * abs(amount)
    return tuple(np.rint(color).astype(np.uint8))


def interpolate_grid(grid: np.ndarray, theta: float, phi: float) -> float:
    theta_step = 180.0 / (grid.shape[0] - 1)
    phi_step = 360.0 / (grid.shape[1] - 1)
    theta_pos = theta / theta_step
    phi_pos = phi / phi_step

    i0 = int(np.floor(theta_pos))
    j0 = int(np.floor(phi_pos))
    i1 = min(i0 + 1, grid.shape[0] - 1)
    j1 = min(j0 + 1, grid.shape[1] - 1)
    theta_frac = theta_pos - i0
    phi_frac = phi_pos - j0

    top = grid[i0, j0] * (1.0 - phi_frac) + grid[i0, j1] * phi_frac
    bottom = grid[i1, j0] * (1.0 - phi_frac) + grid[i1, j1] * phi_frac
    return float(top * (1.0 - theta_frac) + bottom * theta_frac)


def load_font(size: int, bold: bool = False):
    if ImageFont is None:
        return None

    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def draw_centered_text(draw, xy: tuple[int, int], text: str, font, fill=(20, 20, 20)) -> None:
    bbox = draw.textbbox((0, 0), text, font=font)
    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    draw.text((xy[0] - width / 2, xy[1] - height / 2), text, fill=fill, font=font)


def create_difference_map_with_pillow(
    theta_grid: np.ndarray,
    phi_grid: np.ndarray,
    difference_grid: np.ndarray,
    title: str,
    output_path: Path,
) -> None:
    if Image is None:
        raise ModuleNotFoundError("Neither matplotlib nor Pillow is available for PNG output")

    width, height = 2800, 2100
    left, right, top, bottom = 430, 450, 240, 330
    plot_w = width - left - right
    plot_h = height - top - bottom
    colorbar_x = width - right + 115
    colorbar_w = 90

    image = Image.new("RGB", (width, height), "white")
    pixels = image.load()
    vmax = float(np.nanmax(np.abs(difference_grid)))

    for y in range(plot_h):
        theta = 180.0 - (y / (plot_h - 1)) * 180.0
        for x in range(plot_w):
            phi = (x / (plot_w - 1)) * 360.0
            value = interpolate_grid(difference_grid, theta, phi)
            pixels[left + x, top + y] = diverging_color(value, vmax)

    draw = ImageDraw.Draw(image)
    title_font = load_font(76, bold=True)
    label_font = load_font(62)
    tick_font = load_font(52)
    small_font = load_font(46)

    draw_centered_text(draw, (left + plot_w // 2, 105), title, title_font)
    draw_centered_text(draw, (left + plot_w // 2, height - 120), "Azimuthal Angle phi (deg.)", label_font)

    y_label = Image.new("RGBA", (900, 100), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_label)
    y_draw.text((0, 0), "Polar Angle theta (deg.)", fill=(20, 20, 20), font=label_font)
    y_label = y_label.rotate(90, expand=True)
    image.paste(y_label, (55, top + plot_h // 2 - y_label.height // 2), y_label)

    axis_color = (30, 30, 30)
    grid_color = (190, 190, 190)
    for phi_tick in [0, 90, 180, 270, 360]:
        x = left + round((phi_tick / 360.0) * plot_w)
        draw.line([(x, top), (x, top + plot_h)], fill=grid_color, width=2)
        draw.line([(x, top + plot_h), (x, top + plot_h + 18)], fill=axis_color, width=4)
        draw_centered_text(draw, (x, top + plot_h + 70), str(phi_tick), tick_font)

    for theta_tick in [0, 45, 90, 135, 180]:
        y = top + round(((180.0 - theta_tick) / 180.0) * plot_h)
        draw.line([(left, y), (left + plot_w, y)], fill=grid_color, width=2)
        draw.line([(left - 18, y), (left, y)], fill=axis_color, width=4)
        draw_centered_text(draw, (left - 90, y), str(theta_tick), tick_font)

    draw.rectangle([left, top, left + plot_w, top + plot_h], outline=axis_color, width=5)

    for y in range(plot_h):
        value = vmax - (y / (plot_h - 1)) * 2.0 * vmax
        color = diverging_color(value, vmax)
        draw.rectangle([colorbar_x, top + y, colorbar_x + colorbar_w, top + y], fill=color)
    draw.rectangle([colorbar_x, top, colorbar_x + colorbar_w, top + plot_h], outline=axis_color, width=4)

    for fraction, value in [(0.0, vmax), (0.25, vmax / 2.0), (0.5, 0.0), (0.75, -vmax / 2.0), (1.0, -vmax)]:
        y = top + round(fraction * plot_h)
        draw.line([(colorbar_x + colorbar_w, y), (colorbar_x + colorbar_w + 18, y)], fill=axis_color, width=4)
        draw.text((colorbar_x + colorbar_w + 35, y - 25), f"{value:.1f}", fill=axis_color, font=tick_font)

    cbar_label = Image.new("RGBA", (950, 90), (255, 255, 255, 0))
    cbar_draw = ImageDraw.Draw(cbar_label)
    cbar_draw.text((0, 0), "Energy Difference (kcal/mol)", fill=axis_color, font=small_font)
    cbar_label = cbar_label.rotate(90, expand=True)
    image.paste(cbar_label, (width - 105, top + plot_h // 2 - cbar_label.height // 2), cbar_label)

    image.save(output_path)


def create_plots(rows: list[dict[str, float]], output_dir: Path) -> None:
    theta_grid, phi_grid, e_ortho, e_meta, e_para = symmetry_expanded_grids(rows)
    diff_ortho_meta = (e_ortho - e_meta) * AU_TO_KCAL
    diff_para_meta = (e_para - e_meta) * AU_TO_KCAL

    create_difference_map(
        theta_grid,
        phi_grid,
        diff_ortho_meta,
        r"$\Delta E$ (Ortho $-$ Meta)",
        "Delta E (Ortho - Meta)",
        output_dir / "ortho_meta_diff_grid_campaign_no_freq.png",
    )
    create_difference_map(
        theta_grid,
        phi_grid,
        diff_para_meta,
        r"$\Delta E$ (Para $-$ Meta)",
        "Delta E (Para - Meta)",
        output_dir / "para_meta_diff_grid_campaign_no_freq.png",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge grid_campaign_no_freq energies and plot symmetry-expanded maps."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Folder containing the three grid_campaign_no_freq_*_status.csv files.",
    )
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    rows = build_merged_rows(data_dir)
    output_csv = data_dir / OUTPUT_CSV
    write_merged_csv(rows, output_csv)
    create_plots(rows, data_dir)

    print(f"Wrote {len(rows)} merged rows to {output_csv}")
    print(f"Wrote {data_dir / 'ortho_meta_diff_grid_campaign_no_freq.png'}")
    print(f"Wrote {data_dir / 'para_meta_diff_grid_campaign_no_freq.png'}")


if __name__ == "__main__":
    main()
