import argparse
import json
import math
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt

from constants import (
    alpha_lit_dB_cm_MHz,
    c_bone,
    matched_alpha_coeff,
    treeby_actual_attenuation_db_cm,
)


ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
RES = ROOT / "results"
MATLAB_DIR = (Path(__file__).parent / "matlab").as_posix()


def run_matlab_suite(matlab):
    cmd = (
        f"addpath('{MATLAB_DIR}'); "
        f"run_all_benchmarks('{RES.as_posix()}');"
    )
    subprocess.run([matlab, "-batch", cmd], check=True)


def plot_primary_results(data):
    calibration = data["calibration"]
    target = data["constants"]["target_dB_cm_at_f0"]
    values = [
        calibration["uncorrected"]["attenuation_dB_cm"],
        calibration["corrected"]["attenuation_dB_cm"],
        target,
    ]
    errors = [
        calibration["uncorrected"]["relative_error_percent"],
        calibration["corrected"]["relative_error_percent"],
        0.0,
    ]
    positions = [0.0, 0.54, 1.08]
    fig, ax = plt.subplots(figsize=(7, 5.2))
    bars = ax.bar(
        positions,
        values,
        width=0.42,
        color=["#c98b8b", "#7fa8c9", "#d5d0cb"],
        edgecolor="white",
        linewidth=1.0,
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(["uncorrected", "corrected", "target"])
    ax.set_xlim(-0.38, 1.50)
    ax.set_ylabel("Attenuation (dB/cm)")
    ax.set_title("Calibration at 0.5 MHz (interior source, 180 µs)")
    ax.set_ylim(0, max(values) * 1.18)
    ax.set_axisbelow(True)
    ax.grid(axis="y", alpha=0.25)
    for bar, value, error in zip(bars, values, errors):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            f"{value:.4f}\n({error:+.2f}%)",
            ha="center",
            va="bottom",
        )
    fig.tight_layout()
    fig.savefig(FIG / "fig1_calibration_bar.png", dpi=180)
    plt.close(fig)

    sweep = data["frequency_sweep"]
    frequencies = [row["frequency_MHz"] for row in sweep]
    lit_coeff = data["constants"]["alpha_lit_dB_cm_MHz"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(frequencies, [lit_coeff * f for f in frequencies], "--", label="Literature y=1 target")
    ax.plot(frequencies, [row["attenuation_dB_cm"] for row in sweep], "o-", label="alpha_power=2 corrected")
    ax.scatter([0.5], [target], color="crimson", zorder=4)
    corrected_error = calibration["corrected"]["relative_error_percent"]
    ax.annotate(
        f"f0: {corrected_error:+.2f}%",
        (0.5, target),
        xytext=(10, -17),
        textcoords="offset points",
        ha="left",
        va="top",
    )
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Attenuation (dB/cm)")
    ax.set_title("Attenuation vs Frequency: corrected y=2 matches target at f0")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_frequency_sweep.png", dpi=180)
    plt.close(fig)


def plot_exponent_results(all_results, exponents):
    colors = ["#1b9e77", "#377eb8", "#e41a1c", "#984ea3", "#ff7f00"]
    systematic = [0.5, 1.5, 2, 2.1, 2.5]
    diagnostic = [0.9, 1, 1.05, 1.1]
    by_exponent = {round(float(item.get("n", exp)), 4): item for exp, item in zip(exponents, all_results)}

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    for ax, group, title in zip(
        axes,
        [systematic, diagnostic],
        ["Systematic scan", "Near-y=1 diagnostic points (y=1 uses no_dispersion)"],
    ):
        for exponent in group:
            result = by_exponent[round(float(exponent), 4)]
            color = colors[group.index(exponent)]
            sweep = result["frequency_sweep"]
            frequencies = [row["frequency_MHz"] for row in sweep]
            ax.plot(
                frequencies,
                [row["relative_error_percent"] for row in sweep],
                "o-",
                color=color,
                label=f"y={exponent} measured",
            )
            ax.plot(
                frequencies,
                [100 * ((f / 0.5) ** (exponent - 1) - 1) for f in frequencies],
                "--",
                color=color,
                alpha=0.65,
                label=f"y={exponent} ideal power law",
            )
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylabel("Percent error (%)")
        ax.set_title(title)
        ax.grid(alpha=0.25)
        ax.legend(ncol=2, fontsize=8)
    axes[-1].set_xlabel("Frequency (MHz)")
    fig.suptitle("Mismatch depends on exponent, calibration offset, and near-y=1 dispersion")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_error_vs_n.png", dpi=180)
    plt.close(fig)

    lookup = {exp: res for exp, res in by_exponent.items() if not res.get("status")}
    errors_at_1mhz = {
        exp: next(row for row in lookup[exp]["frequency_sweep"] if abs(row["frequency_MHz"] - 1.0) < 1e-6)[
            "relative_error_percent"
        ]
        for exp in lookup
    }
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    x_curve = [0.5 + 2 * i / 299 for i in range(300)]
    ax.plot(x_curve, [100 * ((1 / 0.5) ** (x - 1) - 1) for x in x_curve], "k--", label="ideal power-law ratio")
    ax.plot(systematic, [errors_at_1mhz[round(float(e), 4)] for e in systematic], "o-", color="#377eb8", lw=2, label="systematic scan")
    ax.plot(diagnostic, [errors_at_1mhz[round(float(e), 4)] for e in diagnostic], "D", color="#e66101", ms=7, label="near-y=1 diagnostic")
    ax.margins(x=0.035, y=0.07)
    ax.set_xlabel("Power-law exponent y")
    ax.set_ylabel("Relative error at 1.0 MHz (%)")
    ax.set_title("Exponent summary at 1.0 MHz")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_exponent_summary.png", dpi=180)
    plt.close(fig)


def _analytical_error_curve(y_grid, f0=0.5):
    target = alpha_lit_dB_cm_MHz * f0
    xs, ys = [], []
    for yv in y_grid:
        a0 = matched_alpha_coeff(alpha_lit_dB_cm_MHz, yv, f0)
        att = treeby_actual_attenuation_db_cm(f0, a0, yv, c_bone)
        if math.isfinite(att):
            err = 100.0 * (att - target) / target
            if -95.0 <= err <= 350.0:
                xs.append(yv)
                ys.append(err)
    return xs, ys


def plot_optional_results():
    scan_path = RES / "singularity_scan.json"
    if scan_path.exists():
        scan = sorted(json.loads(scan_path.read_text()), key=lambda row: row["y"])
        stable = [
            row
            for row in scan
            if row.get("field_status", "stable") == "stable"
            and row.get("relative_error_percent") is not None
            and math.isfinite(row["relative_error_percent"])
        ]
        diverged = [row["y"] for row in scan if row.get("field_status") == "diverged"]

        left_clean = [(r["y"], r["relative_error_percent"]) for r in stable if r["y"] < 1.0 and r.get("fit_R2", 1.0) >= 0.99]
        left_sat = [(r["y"], r["relative_error_percent"]) for r in stable if r["y"] < 1.0 and r.get("fit_R2", 1.0) < 0.99]
        right = [(r["y"], r["relative_error_percent"]) for r in stable if r["y"] > 1.0]

        # Exact Treeby & Cox (2014) analytical fractional-Laplacian PDE prediction
        y_ana_left = [0.75 + i * (0.982 - 0.75) / 150 for i in range(151)]
        y_ana_right = [1.006 + i * (1.35 - 1.006) / 150 for i in range(151)]
        ax_l, ay_l = _analytical_error_curve(y_ana_left)
        ax_r, ay_r = _analytical_error_curve(y_ana_right)

        fig, ax = plt.subplots(figsize=(7.8, 5.0))
        ax.axhspan(-5, 5, color="#d9ead3", alpha=0.55, label="within ±5% tolerance")
        if diverged:
            ax.axvspan(min(diverged) - 0.004, 1.0, color="#fddbc7", alpha=0.6,
                       label="Hadamard blow-up (M ≤ 0, |p| > 10³⁰⁰ Pa)")

        ax.plot(ax_l, ay_l, "--", color="#1f78b4", lw=1.6, label="analytical PDE theory (Treeby 2014)")
        ax.plot(ax_r, ay_r, "--", color="#1f78b4", lw=1.6)

        for segment, label in ((left_clean, "k-Wave 1-D measured (R² ≥ 0.999)"), (right, None)):
            if segment:
                ax.plot([p[0] for p in segment], [p[1] for p in segment],
                        "o", color="#b2182b", ms=5.5, zorder=4, label=label)
        if left_sat:
            ax.plot([p[0] for p in left_sat], [p[1] for p in left_sat],
                    "o", mfc="none", mec="#b2182b", mew=1.5, ms=6.5, zorder=4,
                    label="noise-floor saturated (R² < 0.99)")

        ax.axhline(0, color="k", lw=0.8)
        ax.axvline(1.0, color="#555", ls=":", lw=1.2)
        ax.set_xlabel("Power-law exponent y")
        ax.set_ylabel("Calibration error at 0.5 MHz (%)")
        ax.set_title("Near y=1, k-Wave tracks the fractional-Laplacian dispersion pole")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8.2, loc="upper right")
        all_err = [p[1] for p in left_clean + right]
        ax.set_ylim(min(all_err) - 20, max(all_err) + 25)
        ax.set_xlim(0.74, 1.36)
        fig.tight_layout()
        fig.savefig(FIG / "fig5_y1_sensitivity.png", dpi=180)
        plt.close(fig)

    y11_path = RES / "y11_no_dispersion.json"
    y10_path = RES / "y1_no_dispersion_fullband.json"
    if y11_path.exists():
        result = json.loads(y11_path.read_text())
        sweep = result["frequency_sweep"]
        frequencies = [row["frequency_MHz"] for row in sweep]
        simulated = [row["attenuation_dB_cm"] for row in sweep]
        lit_coeff = result["constants"]["alpha_lit_dB_cm_MHz"]
        target = result["constants"]["target_dB_cm_at_f0"]

        fig, ax = plt.subplots(figsize=(7.4, 4.8))
        ax.plot(frequencies, [lit_coeff * f for f in frequencies], "--", color="#555", label="Literature y=1 target")
        if y10_path.exists():
            res10 = json.loads(y10_path.read_text())
            sim10 = [row["attenuation_dB_cm"] for row in res10["frequency_sweep"]]
            ax.plot(frequencies, sim10, "s--", color="#1b9e77", ms=5, lw=1.5, label="y=1.0, no_dispersion (full-band)")
        ax.plot(frequencies, simulated, "o-", color="#377eb8", label="y=1.1, no_dispersion")
        ax.axvline(0.5, color="#555", ls=":", lw=1)
        error_at_f0 = next(row["relative_error_percent"] for row in sweep if abs(row["frequency_MHz"] - 0.5) < 1e-6)
        ax.scatter([0.5], [target], color="crimson", zorder=5)

        ax.set_ylim(min(simulated) * 0.72, max(max(simulated), lit_coeff * frequencies[-1]) * 1.12)
        ax.set_xlim(frequencies[0] - 0.05, frequencies[-1] + 0.04)
        ax.annotate(f"f0 (y=1.1): {error_at_f0:+.2f}%", (0.5, target), xytext=(11, -15),
                    textcoords="offset points", ha="left", va="top", fontsize=9)
        ax.annotate(f"{sweep[0]['relative_error_percent']:+.1f}%",
                    (frequencies[0], simulated[0]), xytext=(12, -4),
                    textcoords="offset points", ha="left", va="top", fontsize=9)
        ax.annotate(f"{sweep[-1]['relative_error_percent']:+.1f}%",
                    (frequencies[-1], simulated[-1]), xytext=(-4, 9),
                    textcoords="offset points", ha="right", va="bottom", fontsize=9)
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Attenuation (dB/cm)")
        ax.set_title("no_dispersion: y=1.0 tracks the band; y=1.1 diverges at band edges")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8.5)
        fig.tight_layout()
        fig.savefig(FIG / "fig6_y11_frequency_scan.png", dpi=180)
        plt.close(fig)


def main(plot_only=False):
    FIG.mkdir(exist_ok=True)
    RES.mkdir(exist_ok=True)
    exponents = [0.5, 0.9, 1, 1.05, 1.1, 1.5, 2, 2.1, 2.5, 3]
    if plot_only:
        plot_primary_results(json.loads((RES / "results.json").read_text()))
        plot_exponent_results(json.loads((RES / "results_n_sweep.json").read_text()), exponents)
        plot_optional_results()
        print(f"[✓] Figures successfully regenerated from pre-computed results in {FIG}.")
        return

    matlab = shutil.which("matlab")
    if not matlab:
        print("[!] MATLAB executable not found on PATH. Falling back to --plot-only mode using archived results.")
        plot_primary_results(json.loads((RES / "results.json").read_text()))
        plot_exponent_results(json.loads((RES / "results_n_sweep.json").read_text()), exponents)
        plot_optional_results()
        print(f"[✓] Figures successfully regenerated from pre-computed results in {FIG}.")
        return

    run_matlab_suite(matlab)
    plot_primary_results(json.loads((RES / "results.json").read_text()))
    plot_exponent_results(json.loads((RES / "results_n_sweep.json").read_text()), exponents)
    plot_optional_results()
    print(f"[✓] Full MATLAB suite and figures completed in {FIG}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run k-Wave experiments or redraw figures from saved JSON results.")
    parser.add_argument("--plot-only", action="store_true", help="redraw figures without launching MATLAB")
    main(plot_only=parser.parse_args().plot_only)
