import json
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def load_result(y):
    tag = str(y).replace(".", "p")
    path = RESULTS / f"near_y_frequency_y{tag}.json"
    if path.exists():
        return json.loads(path.read_text())
    raise FileNotFoundError(path)


def main():
    # Existing complete curves plus the four newly completed near-singularity scans.
    all_results = {}
    # Prefer the committed aggregate for the older curves.
    aggregate = json.loads((RESULTS / "results_n_sweep.json").read_text())
    for item in aggregate:
        if "constants" in item:
            all_results[round(float(item["constants"]["alpha_power"]), 2)] = item
    for y in (0.98, 0.99, 1.01, 1.02):
        all_results[y] = load_result(y)

    colors = {0.90: "#377eb8", 0.98: "#4daf4a", 0.99: "#984ea3", 1.00: "#111111",
              1.01: "#ff7f00", 1.02: "#e41a1c", 1.05: "#a65628", 1.10: "#f781bf"}
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    for y in (0.90, 0.98, 0.99, 1.00, 1.01, 1.02, 1.05, 1.10):
        result = all_results[y]
        sweep = result["frequency_sweep"]
        f = [row["frequency_MHz"] for row in sweep]
        e = [row["relative_error_percent"] for row in sweep]
        style = "s--" if y == 1.00 else "o-"
        label = f"y={y:g}" + (" (no dispersion)" if y == 1.00 else "")
        ax.plot(f, e, style, color=colors[y], lw=2 if y in (0.98, 0.99, 1.00, 1.01, 1.02) else 1.3,
                ms=4, label=label)
    ax.axhspan(-5, 5, color="#d9ead3", alpha=0.45, label="within ±5%")
    ax.axhline(0, color="#444", lw=0.8)
    ax.axvline(0.5, color="#777", ls=":", lw=1)
    ax.set_xlabel("Center frequency (MHz)")
    ax.set_ylabel("Relative attenuation error (%)")
    ax.set_title("Near-y=1 attenuation error across center frequency")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "near_y1_frequency_error.png", dpi=200)
    plt.close(fig)

    # Correct orientation for the original experiment: fixed 0.5 MHz, scan y.
    fixed_frequency = 0.5
    y_values = (0.90, 0.98, 0.99, 1.00, 1.01, 1.02, 1.05, 1.10)
    errors = []
    for y in y_values:
        result = all_results[y]
        row = next(r for r in result["frequency_sweep"]
                   if r["frequency_MHz"] == fixed_frequency)
        errors.append(row["relative_error_percent"])
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.axhspan(-5, 5, color="#d9ead3", alpha=0.45, label="within ±5%")
    ax.plot(y_values, errors, "o-", color="#b2182b", lw=2, ms=6,
            label="fixed frequency = 0.5 MHz")
    ax.axhline(0, color="#444", lw=0.8)
    ax.axvline(1, color="#777", ls=":", lw=1)
    ax.set_xlabel("Power-law exponent y")
    ax.set_ylabel("Relative attenuation error at 0.5 MHz (%)")
    ax.set_title("Fixed 500 kHz: attenuation error versus y")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "fixed_500khz_error_vs_y.png", dpi=200)
    plt.close(fig)

    # Requested experiment view: hold frequency fixed and scan the exponent y.
    scan_frequencies = (0.2, 0.3, 0.4, 0.6, 0.8, 1.0)
    frequency_colors = ("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#e41a1c", "#a65628")
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for frequency, color in zip(scan_frequencies, frequency_colors):
        values = []
        for y in y_values:
            result = all_results[y]
            row = next(r for r in result["frequency_sweep"]
                       if r["frequency_MHz"] == frequency)
            values.append(row["relative_error_percent"])
        for side, indices in enumerate(([0, 1, 2], [4, 5, 6, 7])):
            ax.plot([y_values[i] for i in indices], [values[i] for i in indices],
                    "o-", color=color, lw=1.8, ms=4,
                    label=f"{frequency:g} MHz" if side == 0 else None)
        ax.scatter([1.0], [values[3]], marker="x", color=color, s=35)
    ax.axhspan(-5, 5, color="#d9ead3", alpha=0.45, label="within ±5%")
    ax.axhline(0, color="#444", lw=0.8)
    ax.axvline(1, color="#777", ls=":", lw=1)
    ax.set_xlabel("Power-law exponent y")
    ax.set_ylabel("Relative attenuation error (%)")
    ax.set_title("Fixed-frequency y scans near y=1")
    ax.scatter([], [], marker="x", color="#444", label="y=1: no_dispersion")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "other_frequencies_error_vs_y.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
