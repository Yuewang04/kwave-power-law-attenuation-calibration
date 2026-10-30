import argparse
import json
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
RES = ROOT / "results"
MATLAB_DIR = (Path(__file__).parent / "matlab").as_posix()


def call(matlab, output, exponent):
    subprocess.run(
        [matlab, "-batch", f"addpath('{MATLAB_DIR}'); run_kwave_experiments('{output.as_posix()}',{exponent})"],
        check=True,
    )


def call_controls(matlab, output):
    subprocess.run(
        [matlab, "-batch", f"addpath('{MATLAB_DIR}'); run_controls('{output.as_posix()}')"],
        check=True,
    )


def plot_primary_results(data):
    calibration = data["calibration"]
    target = data["constants"]["target_dB_cm_at_f0"]
    values = [calibration["uncorrected"]["attenuation_dB_cm"], calibration["corrected"]["attenuation_dB_cm"], target]
    errors = [calibration["uncorrected"]["relative_error_percent"], calibration["corrected"]["relative_error_percent"], 0]
    # Bars sit close together and are drawn narrow: the group is compressed on
    # x so the gaps stay small without fattening the bars themselves.
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
    ax.set_title("Calibration at 0.5 MHz")
    ax.set_ylim(0, max(values) * 1.18)
    ax.set_axisbelow(True)
    ax.grid(axis="y", alpha=0.25)
    for bar, value, error in zip(bars, values, errors):
        ax.text(bar.get_x() + bar.get_width() / 2, value, f"{value:.4f}\n({error:+.1f}%)", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(FIG / "fig1_calibration_bar.png", dpi=180)
    plt.close(fig)

    sweep = data["frequency_sweep"]
    frequencies = [row["frequency_MHz"] for row in sweep]
    lit_coeff = data["constants"]["alpha_lit_dB_cm_MHz"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(frequencies, [lit_coeff * f for f in frequencies], "--", label="Literature y=1 target")
    ax.plot(frequencies, [row["attenuation_dB_cm"] for row in sweep], "o-", label="alpha_power=2 corrected")
    ax.scatter([0.5], [target], color="crimson")
    corrected_error = calibration["corrected"]["relative_error_percent"]
    # Below and to the right of f0: both curves pass above that region, so the
    # label does not sit on the dashed target line.
    ax.annotate(f"f0: {corrected_error:+.2f}%", (0.5, target), xytext=(10, -17),
                textcoords="offset points", ha="left", va="top")
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
    by_exponent = {exponent: result for exponent, result in zip(exponents, all_results)}

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    for ax, group, title in zip(axes, [systematic, diagnostic], ["Systematic scan", "Near-y=1 diagnostic points"]):
        for exponent in group:
            result = by_exponent[exponent]
            color = colors[group.index(exponent)]
            sweep = result["frequency_sweep"]
            frequencies = [row["frequency_MHz"] for row in sweep]
            ax.plot(frequencies, [row["relative_error_percent"] for row in sweep], "o-", color=color, label=f"y={exponent} measured")
            ax.plot(frequencies, [100 * ((f / 0.5) ** (exponent - 1) - 1) for f in frequencies], "--", color=color, alpha=0.65, label=f"y={exponent} theory")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylabel("Percent error (%)")
        ax.set_title(title)
        ax.grid(alpha=0.25)
        ax.legend(ncol=2, fontsize=8)
    axes[-1].set_xlabel("Frequency (MHz)")
    fig.suptitle("Mismatch depends on exponent and distance from calibration frequency")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_error_vs_n.png", dpi=180)
    plt.close(fig)

    lookup = {exponent: result for exponent, result in by_exponent.items() if not result.get("status")}
    errors_at_1mhz = {
        exponent: next(row for row in lookup[exponent]["frequency_sweep"] if row["frequency_MHz"] == 1)["relative_error_percent"]
        for exponent in lookup
    }
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    x_curve = [0.5 + 2 * i / 299 for i in range(300)]
    # No ±5% band here: this panel is evaluated at 1.0 MHz, where the theoretical
    # power-law mismatch is itself large away from y=1, so a band derived at the
    # 0.5 MHz calibration frequency would be meaningless. The band belongs to fig5.
    ax.plot(x_curve, [100 * ((1 / 0.5) ** (x - 1) - 1) for x in x_curve], "k--", label="theoretical")
    ax.plot(systematic, [errors_at_1mhz[e] for e in systematic], "o-", color="#377eb8", lw=2, label="systematic scan")
    ax.plot(diagnostic, [errors_at_1mhz[e] for e in diagnostic], "D", color="#e66101", ms=7, label="near-y=1 diagnostic")
    # Margins instead of a hard 0.5-2.5 window: the y=2.5 marker sat half under
    # the right spine when the limits were pinned to the sampled range.
    ax.margins(x=0.035, y=0.07)
    ax.set_xlabel("Power-law exponent y")
    ax.set_ylabel("Relative error at 1.0 MHz (%)")
    ax.set_title("Exponent summary at 1.0 MHz")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_exponent_summary.png", dpi=180)
    plt.close(fig)


def plot_optional_results():
    scan_path = RES / "singularity_scan.json"
    if scan_path.exists():
        scan = sorted(json.loads(scan_path.read_text()), key=lambda row: row["y"])
        x_values = [row["y"] for row in scan]
        y_values = [row["relative_error_percent"] for row in scan]
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        ax.axhspan(-5, 5, color="#d9ead3", alpha=0.55, label="within ±5%")
        # Two segments with a gap at y=1: default dispersion has no valid sample
        # there, and joining across it would imply interpolation over the
        # singularity.
        left = [(x, e) for x, e in zip(x_values, y_values) if x < 1]
        right = [(x, e) for x, e in zip(x_values, y_values) if x > 1]
        for segment, label in ((left, "k-Wave dispersive mode"), (right, None)):
            if segment:
                ax.plot([p[0] for p in segment], [p[1] for p in segment],
                        "o-", color="#b2182b", lw=2, ms=6, label=label)
        ax.axhline(0, color="k", lw=0.8)
        ax.axvline(1, color="#555", ls=":", lw=1)
        ax.set_xlabel("Power-law exponent y")
        ax.set_ylabel("Calibration error at 0.5 MHz (%)")
        ax.set_title("Near y=1, dispersive k-Wave error forms a sharp signed sensitivity zone")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9, loc="best")
        ax.set_ylim(min(y_values) - 15, max(y_values) + 20)
        fig.tight_layout()
        fig.savefig(FIG / "fig5_y1_sensitivity.png", dpi=180)
        plt.close(fig)

    y11_path = RES / "y11_no_dispersion.json"
    if y11_path.exists():
        result = json.loads(y11_path.read_text())
        sweep = result["frequency_sweep"]
        frequencies = [row["frequency_MHz"] for row in sweep]
        simulated = [row["attenuation_dB_cm"] for row in sweep]
        lit_coeff = result["constants"]["alpha_lit_dB_cm_MHz"]
        target = result["constants"]["target_dB_cm_at_f0"]
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        ax.plot(frequencies, [lit_coeff * f for f in frequencies], "--", label="Literature y=1 target")
        ax.plot(frequencies, simulated, "o-", label="y=1.1, no_dispersion")
        ax.axvline(0.5, color="#555", ls=":", lw=1)
        error_at_f0 = next(row["relative_error_percent"] for row in sweep if row["frequency_MHz"] == 0.5)
        ax.scatter([0.5], [target], color="crimson", zorder=3)
        # Labels go into the empty space beside each point rather than on top of
        # it: the two curves converge at f0, so any label placed at the point
        # itself collides with a line.
        ax.set_ylim(min(simulated) * 0.72, max(max(simulated), lit_coeff * frequencies[-1]) * 1.12)
        ax.set_xlim(frequencies[0] - 0.05, frequencies[-1] + 0.04)
        ax.annotate(f"f0: {error_at_f0:+.2f}%", (0.5, target), xytext=(11, -15),
                    textcoords="offset points", ha="left", va="top")
        ax.annotate(f"{sweep[0]['relative_error_percent']:+.1f}%",
                    (frequencies[0], simulated[0]), xytext=(12, 0),
                    textcoords="offset points", ha="left", va="center")
        ax.annotate(f"{sweep[-1]['relative_error_percent']:+.1f}%",
                    (frequencies[-1], simulated[-1]), xytext=(-3, 9),
                    textcoords="offset points", ha="right", va="bottom")
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Attenuation (dB/cm)")
        ax.set_title("y=1.1 no_dispersion: matched at f0, not across the full band")
        ax.grid(alpha=0.25)
        ax.legend()
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
        return

    matlab = shutil.which("matlab")
    if not matlab:
        raise RuntimeError("MATLAB executable not found on PATH; add MATLAB to PATH before running.")
    raw = RES / "matlab_results.json"
    call(matlab, raw, 2)
    data = json.loads(raw.read_text())
    (RES / "results.json").write_text(json.dumps(data, indent=2))
    plot_primary_results(data)

    all_results = []
    for exponent in exponents:
        output = RES / f"_n_{str(exponent).replace('.', 'p')}.json"
        try:
            call(matlab, output, exponent)
            result = json.loads(output.read_text())
            result["n"] = exponent
            all_results.append(result)
            output.unlink()
        except subprocess.CalledProcessError:
            # The exact rejection message is not captured by this batch driver.
            # alpha_power=3 is an odd integer and tan(pi*y/2) is singular at every
            # odd integer, so y=3 shares the y=1 singularity class. Run the MATLAB
            # call by hand to see the message k-Wave actually returns.
            all_results.append({"n": exponent, "status": "unsupported_by_kwave",
                "reason": "MATLAB/k-Wave rejected this alpha_power; the exact message is not captured by the batch driver. Note tan(pi*y/2) is singular at every odd integer, so y=3 shares the y=1 singularity class."})
    (RES / "results_n_sweep.json").write_text(json.dumps(all_results, indent=2))
    plot_exponent_results(all_results, exponents)
    plot_optional_results()
    call_controls(matlab, RES / "controls_y098.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run k-Wave experiments or redraw figures from saved JSON results.")
    parser.add_argument("--plot-only", action="store_true", help="redraw figures without launching MATLAB")
    main(plot_only=parser.parse_args().plot_only)
