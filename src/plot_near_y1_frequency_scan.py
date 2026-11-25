import json
import math
from pathlib import Path

import matplotlib.pyplot as plt

from constants import (
    alpha_lit_dB_cm_MHz,
    c_bone,
    matched_alpha_coeff,
    treeby_actual_attenuation_db_cm,
)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def load_result(y):
    tag = str(y).replace(".", "p")
    path = RESULTS / f"near_y_frequency_y{tag}.json"
    if path.exists():
        return json.loads(path.read_text())
    raise FileNotFoundError(path)


def _is_valid_row(row):
    val = row.get("relative_error_percent")
    return (
        row.get("field_status", "stable") == "stable"
        and val is not None
        and math.isfinite(val)
    )


def main():
    all_results = {}
    aggregate = json.loads((RESULTS / "results_n_sweep.json").read_text())
    for item in aggregate:
        if "constants" in item:
            all_results[round(float(item["constants"]["alpha_power"]), 2)] = item
    for y in (0.98, 0.99, 1.01, 1.02):
        all_results[y] = load_result(y)

    # 1. Frequency error curves for STABLE near-y=1 exponents (y=0.99 excluded due to Hadamard blow-up)
    colors = {
        0.90: "#377eb8",
        0.98: "#4daf4a",
        1.00: "#111111",
        1.01: "#ff7f00",
        1.02: "#e41a1c",
        1.05: "#a65628",
        1.10: "#f781bf",
    }
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    for y in (0.90, 0.98, 1.00, 1.01, 1.02, 1.05, 1.10):
        result = all_results[y]
        sweep = [r for r in result["frequency_sweep"] if _is_valid_row(r)]
        f = [row["frequency_MHz"] for row in sweep]
        e = [row["relative_error_percent"] for row in sweep]
        style = "s--" if y == 1.00 else "o-"
        label = f"y={y:g}" + (" (no_dispersion)" if y == 1.00 else "")
        ax.plot(
            f,
            e,
            style,
            color=colors[y],
            lw=2.0 if y in (0.98, 1.00, 1.01, 1.02) else 1.4,
            ms=4.5,
            label=label,
        )
    ax.axhspan(-5, 5, color="#d9ead3", alpha=0.45, label="within ±5%")
    ax.axhline(0, color="#444", lw=0.8)
    ax.axvline(0.5, color="#777", ls=":", lw=1)
    ax.set_xlabel("Center frequency (MHz)")
    ax.set_ylabel("Relative attenuation error (%)")
    ax.set_title("Near-y=1 stable exponents across frequency (y=0.99 excluded: Hadamard blow-up)")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "near_y1_frequency_error.png", dpi=200)
    plt.close(fig)

    # 2. Fixed 500 kHz scan across y: NEVER connect across y=1 or through diverged y=0.99!
    fixed_frequency = 0.5
    left_y = (0.90, 0.98)
    right_y = (1.01, 1.02, 1.05, 1.10)

    def get_err(yv, freq):
        row = next(r for r in all_results[yv]["frequency_sweep"] if abs(r["frequency_MHz"] - freq) < 1e-6)
        return row["relative_error_percent"] if _is_valid_row(row) else None

    target_f0 = alpha_lit_dB_cm_MHz * fixed_frequency
    ana_l_y = [0.90 + i * (0.984 - 0.90) / 80 for i in range(81)]
    ana_r_y = [1.008 + i * (1.10 - 1.008) / 80 for i in range(81)]
    ana_l_e = [100.0 * (treeby_actual_attenuation_db_cm(fixed_frequency, matched_alpha_coeff(alpha_lit_dB_cm_MHz, yv, fixed_frequency), yv, c_bone) - target_f0) / target_f0 for yv in ana_l_y]
    ana_r_e = [100.0 * (treeby_actual_attenuation_db_cm(fixed_frequency, matched_alpha_coeff(alpha_lit_dB_cm_MHz, yv, fixed_frequency), yv, c_bone) - target_f0) / target_f0 for yv in ana_r_y]

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.axhspan(-5, 5, color="#d9ead3", alpha=0.45, label="within ±5%")
    ax.axvspan(0.986, 1.00, color="#fddbc7", alpha=0.65, label="y=0.99 blow-up (M ≤ 0)")
    ax.plot(ana_l_y, ana_l_e, "--", color="#1f78b4", lw=1.6, label="analytical PDE theory")
    ax.plot(ana_r_y, ana_r_e, "--", color="#1f78b4", lw=1.6)
    ax.plot(left_y, [get_err(y, fixed_frequency) for y in left_y], "o", color="#b2182b", ms=6.5, label="k-Wave default (stable)")
    ax.plot(right_y, [get_err(y, fixed_frequency) for y in right_y], "o", color="#b2182b", ms=6.5)
    err_100 = get_err(1.00, fixed_frequency)
    ax.scatter([1.00], [err_100], marker="x", color="#111111", s=65, lw=2, zorder=5, label="y=1.00 (no_dispersion)")
    ax.axhline(0, color="#444", lw=0.8)
    ax.axvline(1.0, color="#777", ls=":", lw=1)
    ax.set_xlabel("Power-law exponent y")
    ax.set_ylabel("Relative attenuation error at 0.5 MHz (%)")
    ax.set_title("Fixed 500 kHz: default dispersion vs analytical PDE theory (no cross-singularity join)")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIGURES / "fixed_500khz_error_vs_y.png", dpi=200)
    plt.close(fig)

    # 3. Other fixed frequencies: plot ONLY stable points, split left/right of y=1, mark y=1.0 as 'x'
    scan_frequencies = (0.2, 0.3, 0.4, 0.6, 0.8, 1.0)
    frequency_colors = ("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#e41a1c", "#a65628")
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    ax.axvspan(0.986, 1.00, color="#fddbc7", alpha=0.55, label="y=0.99 blow-up (excluded)")
    for frequency, color in zip(scan_frequencies, frequency_colors):
        left_pts = [(y, get_err(y, frequency)) for y in left_y if get_err(y, frequency) is not None]
        right_pts = [(y, get_err(y, frequency)) for y in right_y if get_err(y, frequency) is not None]
        ax.plot([p[0] for p in left_pts], [p[1] for p in left_pts], "o-", color=color, lw=1.8, ms=4.5, label=f"{frequency:g} MHz")
        ax.plot([p[0] for p in right_pts], [p[1] for p in right_pts], "o-", color=color, lw=1.8, ms=4.5)
        e1 = get_err(1.00, frequency)
        if e1 is not None:
            ax.scatter([1.0], [e1], marker="x", color=color, s=40, zorder=5)
    ax.axhspan(-5, 5, color="#d9ead3", alpha=0.45, label="within ±5%")
    ax.axhline(0, color="#444", lw=0.8)
    ax.axvline(1, color="#777", ls=":", lw=1)
    ax.set_xlabel("Power-law exponent y")
    ax.set_ylabel("Relative attenuation error (%)")
    ax.set_title("Fixed-frequency y scans near y=1 (diverged y=0.99 excluded)")
    ax.scatter([], [], marker="x", color="#444", label="y=1: no_dispersion")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "other_frequencies_error_vs_y.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
