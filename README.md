# kwave-power-law-attenuation-calibration

[中文 README](README_ZH.md)

[Experimental report (PDF)](report.pdf). LaTeX source: [`report.tex`](report.tex); rebuild with `pdflatex report.tex`, three passes.

> **Key idea:** Correctly converting a literature attenuation coefficient to k-Wave's power-law parameterization does not by itself guarantee correct numerical attenuation. This project separates those two questions and audits both under a controlled 1-D setup.

![MATLAB R2026a](https://img.shields.io/badge/MATLAB-R2026a-e86e05) ![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-3776ab) [![License: MIT](https://img.shields.io/badge/License-MIT-2ea44f)](LICENSE)

## Contents

- [TL;DR](#tldr)
- [Background and motivation](#background-and-motivation)
- [Quick guidance: which y should I use?](#quick-guidance-which-y-should-i-use)
- [General formula and derivation](#general-formula-and-derivation)
- [Part 1: Analytical calibration and exponent sweep](#part-1-analytical-calibration-and-exponent-sweep)
- [Part 2: Near-y=1 sensitivity of peak-pressure attenuation estimates](#part-2-near-y1-sensitivity-of-peak-pressure-attenuation-estimates)
- [Relation to existing k-Wave guidance](#relation-to-existing-k-wave-guidance)
- [Part 3: Practical implication and worked example](#part-3-practical-implication-and-worked-example)
- [Conclusion](#conclusion)
- [Reproduction](#reproduction)
- [Limitations](#limitations)
- [Author](#author)
- [Citation](#citation)
- [References](#references)
- [Feedback and license](#feedback-and-license)

## TL;DR

- **Problem 1:** If a literature coefficient reported for a linear frequency law (y=1) is used directly in k-Wave's y=2 power-law model without adjustment, the simulated attenuation is biased away from the calibration frequency—in this project's y=2, 0.5-MHz test, by −50.13%. Matching the coefficient with the derived formula reduced that to −0.97%.
- **Problem 2:** Using a power-law exponent close to but not equal to 1 under k-Wave's default dispersion mode can silently produce large calibration errors, even though the simulation runs without error and reports a high R². In this project's tests, that error exceeded ±5% across a broad y≈1 region (sampled exponents outside the band ran from y=0.80, +5.48%, to y=1.32, −5.07%); switching to `no_dispersion` at y=1.1 reduced the error to −0.50%.
- **The near-y=1 error is not a power-law mismatch.** At y=1.01 the attenuation deficit is −56.8% to −59.4% across 0.2–1.0 MHz with R²≥0.999 — nearly independent of frequency, whereas the theoretical power-law mismatch over that band is at most ±1.4%. A frequency-independent offset is a multiplicative attenuation deficit in the solver, not a wrong exponent. An independent 3-D control shows the matching signature: the same y=1.01 configuration raises peak pressure by +57% (2.74 vs 1.75 MPa) relative to `no_dispersion`.
- **Near y=1 the field itself can fail, and a percent-error table hides it.** At y=0.99 the 1-D peak pressure collapses to the 10⁻⁵ Pa level and rises with distance, so no attenuation estimate exists for that exponent; in the 3-D control the same configuration diverges to ~10³⁰ Pa in both the MATLAB and CUDA paths at every timestep tested. These are qualitative failures, not calibration biases.

## Background and motivation

Transcranial ultrasound simulations need reliable skull attenuation coefficients. Two separate questions arise: literature values are often reported with a linear frequency law (y=1), whereas k-Wave uses a configurable power law, so a coefficient copied directly into a different exponent is not matched away from the reference frequency; and the dispersion relation underlying k-Wave's power-law absorption model is singular at y=1, so accepted values near 1 can produce substantial deviations between intended and simulated attenuation.

I came across this problem while preparing an acoustic simulation dataset. The benchmark I was following used y=2, while the bone attenuation values I found in the literature used y=1. I first treated that as a parameter-definition detail, then built a small 1-D benchmark when the conversion did not behave as expected.

In the project's 1-D test, an uncorrected coefficient biased the simulated pressure attenuation by roughly 50%. That magnitude is worth checking independently before using k-Wave output for dosing or focal planning. The benchmark is controlled 1-D evidence, not a result for real skull geometry.

## Quick guidance: which y should I use?

Use the route that matches your simulation goal:

**Route A — y=2, default dispersion.** This is the convention used by the Aubry et al. (2022) benchmark and the combination fully validated in Part 1. It suits most workflows, especially when retaining frequency-dependent phase coupling matters (as it does more for broadband or pulsed simulations). Use `alpha_coeff = lit_coeff × f0^(1-y)` with y=2. **f0 must be entered in MHz, not Hz**; because the exponent is generally non-integer, the result is not a simple power-of-ten conversion. In k-Wave's power-law absorption formulation, the dispersion correction contains tan(πy/2), which is zero at y=2. Thus `alpha_mode` does not change that specific term at y=2; this is not a general guarantee that every solver or backend is equivalent.

**Route B — y=1 with `no_dispersion` (validated here only with the MATLAB solver).** k-Wave has a mathematical singularity at y=1, so default dispersion reports an error and `alpha_mode='no_dispersion'` is required. The derived ratio `(f/f0)^(y-1)` is identically 1 at y=1, so this combination should theoretically match the literature target over the full band. This project validated the calibration point (−0.45%), but did not run a dedicated full-band scan as it did for y=2 or y=1.1. Applications that require accurate phase calculations, including transcranial focusing plans and phase-aberration correction, should assess Route B carefully because removing frequency-dependent sound-speed coupling may affect focusing accuracy. ⚠️ If you plan to use a GPU/CUDA backend, see the backend portability warning in [Relation to existing k-Wave guidance](#relation-to-existing-k-wave-guidance) before relying on `no_dispersion`.

**Route C — y in 0.80–1.32, excluding 1.** Part 2 found default-dispersion calibration errors ranging from the tens of percent to over 100% in this interval, even when the program accepts the input and R² remains high. Unless you have read Part 2 and understand this trap, avoid default dispersion in this range.

**No Route D is recommended near y=1.** Avoiding exact equality with 1 does not by itself establish a safe workaround. The repository does not contain a standalone saved result file for the y=0.99 timestep-refinement and y=1.01 attenuation-estimator checks, so those specific follow-up numbers are not reported here. No universal safe distance from y=1 has been established.

### Before/after code: matching the coefficient

The diff below replaces the common uncorrected assignments with the reference-frequency conversion. `f0_MHz` is explicitly in **MHz**.

```diff
 f0_MHz = 0.5;                 % MHz (do not use Hz here)
 lit_coeff = 1.86749;          % dB/(cm*MHz), y=1 literature law
-medium.alpha_coeff = lit_coeff;                    % WRONG: copied directly
-medium.alpha_power = 2;
+y = 2;
+medium.alpha_coeff = lit_coeff * f0_MHz^(1-y);     % RIGHT: match at f0
+medium.alpha_power = y;
```

## General formula and derivation

Let the literature law be $\alpha_{\mathrm{lit}}(f)=\alpha_{\mathrm{lit,coeff}}f$, and k-Wave's law be $\alpha_{\mathrm{kwave}}(f)=\alpha_{\mathrm{coeff}}(y)f^y$. Matching at calibration frequency $f_0$ requires

$$
\alpha_{\mathrm{coeff}}(y)f_0^y=\alpha_{\mathrm{lit,coeff}}f_0^1.
$$

Dividing by f0^y gives

$$
\alpha_{\mathrm{coeff}}(y)=\alpha_{\mathrm{lit,coeff}}f_0^{1-y},\qquad
\frac{\alpha_{\mathrm{kwave}}(f)}{\alpha_{\mathrm{lit}}(f)}=\left(\frac{f}{f_0}\right)^{y-1}.
$$

Note: f0 must be expressed in **MHz** to match k-Wave's native `alpha_coeff` unit convention, dB/(cm·MHz^y). Substituting f0 in **Hz** will not simply rescale by a power of ten: the exponent (1-y) is generally non-integer.

At y=1, the exponent in the ratio is zero, so the ratio is 1 at every frequency. Thus y=1 with `no_dispersion` should theoretically track the literature law across the full band, although this project validated only the calibration point (−0.45% error, R²=0.999994), not a dedicated full-band y=1 scan. Exactly y=1 is a special case for k-Wave's default dispersion mode and requires `alpha_mode='no_dispersion'`. The y=0.9, 1.05, and 1.1 cases are diagnostic samples for the near-y=1 sensitivity, not recommended operating values. Default dispersion retains causal Kramers–Kronig phase coupling, which is why benchmark workflows often use y=2; y=1 with `no_dispersion` gives up that frequency-dependent sound-speed coupling.

## Part 1: Analytical calibration and exponent sweep

The 1-D MATLAB k-Wave experiment used `alpha_lit_coeff=1.86749 dB/(cm·MHz)`, f0=0.5 MHz, c=3100 m/s, rho=2200 kg/m³, an 18-cycle toneburst, and six sensors spanning 26–126 mm. At each sensor, the recorded `p_max` is the maximum positive pressure returned by k-Wave; the code takes its absolute value before fitting ln(p_max) versus distance. Thus the reported quantity is a peak-pressure-based attenuation estimate, not a direct measurement of the material's intrinsic attenuation coefficient. The systematic points were y=0.5, 1.5, 2, 2.1, and 2.5. We also sampled y=0.9, 1.05, and 1.1 as diagnostic points near y=1.
All six source-to-sensor distances are included in the fit. The y=1 row uses `no_dispersion` as a special-case reference and is not a like-for-like default-dispersion comparison.

Read the table: `calibration error` is the percent difference between the simulated attenuation and the literature target at f0. `status` marks `usable (sensitivity)` for the three y values sampled specifically to probe the near-y=1 instability described in Part 2—these are not recommended operating exponents, only diagnostic points.

| y | alpha_coeff | calibration error (%) | R² | status |
|---:|---:|---:|---:|---|
| 0.5 | 1.320515 | +1.06 | 0.999995 | usable |
| 0.9 | 1.742430 | +13.59 | 0.999999 | usable (sensitivity) |
| 1.0 | 1.867490 | −0.45 | 0.999994 | usable (`no_dispersion`) |
| 1.05 | 1.933347 | −22.12 | 0.999992 | usable (sensitivity) |
| 1.1 | 2.001526 | −13.15 | 0.999999 | usable (sensitivity) |
| 1.5 | 2.641030 | −3.28 | 0.999996 | usable |
| 2.0 | 3.734980 | −0.97 | 0.999986 | usable |
| 2.1 | 4.003052 | −0.38 | 0.999989 | usable |
| 2.5 | 5.282059 | +3.73 | 0.999985 | usable |
| 3.0 | — | — | — | `unsupported_by_kwave` |

y=3.0 was rejected by k-Wave and no attenuation value is reported for it. Note that the dispersion factor `tan(πy/2)` is singular at every odd integer, so y=3 belongs to the same singularity class as y=1; the rejection message itself was not preserved in the committed results, so the exact reason k-Wave gave is not asserted here.

<img src="figures/fig1_calibration_bar.png" width="560" alt="Bar chart comparing uncorrected, corrected, and target attenuation at 0.5 MHz">

Figure 1. Bar chart at 0.5 MHz. Three bars show attenuation in dB/cm: Uncorrected (literature coefficient used directly with y=2), Corrected (coefficient matched to f0 using this project's formula), and Target (the literature y=1 value at f0). Each bar is labeled with its value and percent error versus Target. Matching removes the large calibration bias seen in the Uncorrected bar.

<img src="figures/fig2_frequency_sweep.png" width="560" alt="Attenuation versus frequency showing corrected y=2 model and literature target">

Figure 2. Line chart. X-axis: frequency (MHz), 0.2 to 1.0. Y-axis: attenuation (dB/cm). Dashed line: the literature y=1 target at each frequency. Solid line: the corrected y=2 model's simulated attenuation at each frequency. The two lines meet only at f0=0.5 MHz (marked); away from f0, the y=2 model's attenuation grows faster with frequency than the literature's y=1 line, which is the expected consequence of the exponents differing.

<img src="figures/fig3_error_vs_n.png" width="560" alt="Measured and theoretical relative attenuation errors across frequencies and exponents">

Figure 3. Two-panel line chart. X-axis (both panels): frequency (MHz). Y-axis (both panels): percent error versus the literature target. Each exponent y has two lines of the same color: solid = measured (simulated) error, dashed = theoretical error predicted by (f/f0)^(y-1)-1. Top panel: systematic scan points (y=0.5, 1.5, 2, 2.1, 2.5), where measured and theoretical lines closely overlap. Bottom panel: near-y=1 diagnostic points (y=0.9, 1, 1.05, 1.1), where the measured lines for y=1.05 and y=1.1 visibly depart from their theoretical dashed counterparts.

<img src="figures/fig4_exponent_summary.png" width="560" alt="Exponent summary showing systematic trend and near-y=1 sensitivity dip">

Figure 4. Scatter/line chart. X-axis: power-law exponent y, 0.5 to 2.5. Y-axis: relative error (%) at 1.0 MHz. Black dashed line: ideal theoretical error curve. Blue connected dots: systematic scan points, which closely follow the theoretical curve. Orange diamonds (not connected by a line): near-y=1 diagnostic points, which fall well off the theoretical curve. At 1.0 MHz the theoretical curve is itself far from zero away from y=1, so this panel shows departure from the expected power-law mismatch rather than a ±5% calibration band; that band is shown only in Figure 5, where it is evaluated at the calibration frequency.

The algebraic conversion answers a parameter-transfer question. It does not guarantee that a numerical solver will reproduce the converted attenuation. One practical warning is worth keeping in view: a high R² only indicates approximately exponential pressure decay over the fitted interval. It does not certify the recovered attenuation coefficient.

<a id="part-2-near-y1-sensitivity-of-peak-pressure-attenuation-estimates"></a>
## Part 2: Near-y=1 sensitivity of peak-pressure attenuation estimates

k-Wave documents y=1 as a special input case and suggests a nearby exponent or `alpha_mode='no_dispersion'`. This mode disables the dispersion correction that couples power-law absorption to causal Kramers–Kronig phase behavior. That is reasonable for this narrowband toneburst, but it gives up frequency-dependent sound-speed coupling and can distort phase in broadband or pulse simulations. We first encountered the issue as an unexpected 0.2-MHz result. After checking the main boundary-related settings, the same pattern appeared while extending the test across y values. The fine scan then reached +143.11% at y=0.98 (R²=0.96) and −41.01% at y=1.02. The sampled exponents outside ±5% ran from y=0.80 (+5.48%) to y=1.32 (−5.07%); the adjacent in-band samples were y=0.78 (+4.83%) and y=1.35 (−4.71%). These are setup-dependent reference boundaries, not universal constants. The reported quantity is a peak-pressure attenuation estimate; the present tests do not fully separate dispersion/model effects, waveform reshaping, estimator bias, and numerical discretization.

### Control checks

For the representative y=0.98 case, timestep refinement gave peak-pressure relative errors of 143.11% (baseline), 143.72% (dt/2), and 143.79% (dt/4), with R² values of 0.9609, 0.9614, and 0.9616. The anomaly changed little under this refinement test. A frequency-domain amplitude fit around the nearest f0 bin gave 210.10%, 181.49%, and 181.40% for the same three runs (R²=0.5797, 0.6659, and 0.6663). The difference shows that the sensitivity depends on the attenuation estimator; these controls do not identify a unique intrinsic material attenuation.

These control values are stored in [`results/controls_y098.json`](results/controls_y098.json), generated by [`src/matlab/run_controls.m`](src/matlab/run_controls.m) through the Python batch driver.

<img src="figures/fig5_y1_sensitivity.png" width="560" alt="Calibration error versus exponent showing sensitivity near y=1">

Figure 5. Line chart, fine scan. X-axis: power-law exponent y, 0.75 to 1.35. Y-axis: calibration error (%) at 0.5 MHz using default dispersion mode. Shaded green band: within ±5% of target. Vertical dotted line: y=1. The error spikes sharply and changes sign on either side of y=1. The curve is drawn as two segments that are deliberately not joined across y=1, because default dispersion has no valid sample there; no interpolation across the singularity is implied.

y=1 is a special case, so its reference point uses `no_dispersion` and is excluded from like-for-like comparisons with default-dispersion points.

### Amplitude dependence of the attenuation estimate

A calibration error is normally read as a fixed multiplicative bias. Near y=1 it is not one. Holding y=1.05 fixed and scaling `alpha_coeff` by 0.5×, 1×, and 2× gives measured attenuations of 0.407, 0.727, and 1.189 dB/cm against power-law expectations of 0.451, 0.902, and 1.804 dB/cm ([data](results/tan_scaling_check.json)):

| `alpha_coeff` scale | `alpha_coeff` | Measured (dB/cm) | Power-law expectation (dB/cm) | Error |
|---|---:|---:|---:|---:|
| 0.5× | 0.9667 | 0.407 | 0.451 | −9.8% |
| 1× | 1.9333 | 0.727 | 0.902 | −19.4% |
| 2× | 3.8667 | 1.189 | 1.804 | −34.1% |

Doubling the coefficient multiplies the fitted attenuation by 1.79× and then 1.64×, not by 2×. For a wave decaying as exp(−αx) the fitted slope must be proportional to `alpha_coeff` whatever the source amplitude, so this compression means the recorded pressure is not decaying exponentially over the fitted interval. Grid refinement does not remove it: at the same exponent, dx = 0.5, 0.25, and 0.125 mm give −22.12%, −21.86%, and −21.72% ([data](results/grid_resolution_check.json)), so it is not a discretization artifact.

Two consequences follow. The percentages quoted in this section are conditional on the coefficient used — "−22.12% at y=1.05" is the value at `alpha_coeff = 1.9333`, not a property of the exponent alone. And the near-y=1 estimate is better read as a distorted peak-pressure observable than as a material attenuation coefficient, which is consistent with the collapse at y=0.99. The same scaling check has not been run at other exponents, so whether this compression is specific to y≈1 or also affects the systematic scan is not established here; that is the most useful follow-up experiment.

### Fixed-frequency y scans beyond 500 kHz

The additional scans hold the center frequency at 0.2, 0.3, 0.4, 0.6, 0.8, or 1.0 MHz and vary y. New runs add y=0.98, 0.99, 1.01, and 1.02 to the existing y=0.9, 1.0, 1.05, and 1.1 results. All coefficients remain matched at f0=0.5 MHz; the reference is the linear literature attenuation at each tested frequency. Away from f0, the plotted error therefore includes the expected power-law mismatch, `(f/f0)^(y-1)-1`, as well as deviations of the peak-pressure estimate.

<img src="figures/other_frequencies_error_vs_y.png" width="760" alt="Fixed-frequency attenuation-error curves versus y at 0.2, 0.3, 0.4, 0.6, 0.8, and 1.0 MHz; y=1 is a separate no-dispersion reference">

Supplementary Figure S1. Each curve holds frequency fixed and scans y. Connected circles use default dispersion, including y=1.1; separate crosses at y=1 use `no_dispersion`. Lines connect sampled points for visual guidance and do not interpolate across y=1. The shaded band marks ±5% relative to the linear literature target.

The observed sensitivity extends beyond 500 kHz. At y=0.98, errors range from +485.2% at 0.2 MHz to −7.2% at 1.0 MHz. At y=0.99 the recorded peak pressure collapsed to the 10⁻⁵ Pa level (1.7×10⁻⁵ to 7.7×10⁻⁵ Pa across the six frequencies) and increased with distance (at 0.5 MHz, 7.3×10⁻⁵ to 9.8×10⁻⁵ Pa from 26 to 126 mm, R²=0.69–0.83), so the fitted slope describes numerical residue rather than a decaying wave; no attenuation estimate is reported for this exponent. At y=1.01 and 1.02, errors are approximately −57% to −59% and −41% to −42%, and are nearly flat across the 0.2–1.0 MHz range, where the theoretical power-law mismatch is at most ±1.4%. A frequency-independent offset of that size is a multiplicative attenuation deficit rather than a mismatch in the frequency exponent.

These are exploratory results from the original setup — the same configuration as Part 1, in which the source sits inside the PML — with a 90-µs record and dx=0.5 mm. The 0.2-MHz burst lasts approximately 90 µs, creating a record-truncation risk, and some near-y=1 runs emit a timestep warning. Convergence and estimator controls have not been completed across this matrix. The results do not establish a universal failure interval or a causal explanation.

Data: [`y=0.98`](results/near_y_frequency_y0p98.json), [`y=0.99`](results/near_y_frequency_y0p99.json), [`y=1.01`](results/near_y_frequency_y1p01.json), [`y=1.02`](results/near_y_frequency_y1p02.json), and the existing [`exponent results`](results/results_n_sweep.json). Run `run_near_y_frequency_scans` from [`src/matlab/run_near_y_frequency_scans.m`](src/matlab/run_near_y_frequency_scans.m), then `python src/plot_near_y1_frequency_scan.py` to regenerate the figure. That script also writes two auxiliary views, `figures/near_y1_frequency_error.png` (error versus frequency) and `figures/fixed_500khz_error_vs_y.png`, which are committed alongside it. See the [experiment design](experimental_design_near_y1.md) for controls still required.

## Relation to existing k-Wave guidance

k-Wave documentation already identifies y=1 as a special case, provides `no_dispersion` as an option, and discusses reference-frequency adjustment of `alpha_coeff`. Here, those pieces are tested together under one controlled 1-D setup. This work does not claim novelty for these individual pieces.

Informal k-Wave User Forum discussions provide useful historical context, but they are not peer-reviewed evidence and their version-specific applicability cannot be confirmed here:

- In a thread dated roughly 2016, one user reported that adding `medium.alpha_mode = 'no_dispersion'` made no difference to that user's simulation results ([thread](https://www.k-wave.org/forum/topic/attenuation-in-heroginious-medium)).
- In a later thread, a user explicitly working with k-Wave 1.2 reported that the C++/CUDA path ignored `no_dispersion`, while the MATLAB path used it; the same report described no 3-D difference but a difference in 1-D and suspected an old-version bug ([thread](https://www.k-wave.org/forum/topic/attenuation-measurement-in-3d)). A separate developer reply also stated that the C++/CUDA code did not support `alpha_mode` ([thread](https://www.k-wave.org/forum/topic/vastly-different-time-of-travel-for-the-same-model-in-3d-and-2d)).
- In another discussion, a k-Wave developer explained that `no_dispersion` removes one term in the absorption operator: its effect is strongest when the exponent is close to 1 and weaker farther away ([thread](https://www.k-wave.org/forum/topic/relationship-between-attenuation-and-wave-velocity)).

The second report concerns k-Wave 1.2. This repository's core 1-D benchmark uses the MATLAB solver. A separate small 3-D backend control found that the tested k-Wave 1.4.1 wrapper/CUDA v1.3 pairing did not apply `no_dispersion` and that the CUDA binary rejected exact y=1 ([result](results/cuda_y1_backend_control_20260914.json)). In the same 3-D control, y=0.99 under default dispersion produced peak pressures of order 10³⁰ Pa in both the MATLAB and the CUDA path, and stayed of that order at CFL 0.1, 0.05, and 0.025 — a nonphysical field rather than an attenuation value, and one that timestep refinement does not remove. This is a version-specific control, not a claim about every CUDA release. Readers should test `no_dispersion` in their own backend, dimensionality, and toolbox version.

Backend portability warning: verifying `no_dispersion` in the MATLAB solver does not establish that a GPU/CUDA backend receives or applies the same setting. Before changing an exponent away from a well-tested value, compare default and `no_dispersion` outputs on the exact backend used for production.

## Part 3: Practical implication and worked example

### Why this matters in a CT-based pipeline

The calibration issue becomes practical when a CT-to-acoustic-property mapper computes an attenuation prefactor and passes it to k-Wave. The function below is written from scratch to demonstrate this class of error—it isn't lifted from, or meant to represent, any particular codebase.

```matlab
% WARNING: illustrative uncorrected mistake pattern; do not copy to production.
function medium = map_hu_to_medium_properties(HUmat)
    psi = max(0, min(1, (1500 - single(HUmat)) / 1500));
    medium.sound_speed = 1483*psi + 2800*(1-psi);
    medium.density = 1000*psi + 1900*(1-psi);
    % WARNING: synthetic HU-to-property mapping; not a validated CT model.
    synthetic_alpha_prefactor_np_per_m_mhz = 21.5 + ...
        (208.9 - 21.5) .* sqrt(psi);
    medium.alpha_coeff = synthetic_alpha_prefactor_np_per_m_mhz * 0.08686;
    medium.alpha_power = 1.1;
    % alpha_mode='no_dispersion' is omitted in the mistake pattern
end
```

The function above represents a complete synthetic CT-to-medium mapping pattern: sound speed, density, and attenuation vary spatially with `psi`. If it were used directly in a heterogeneous simulation, spatial changes in sound speed and density would produce impedance-mismatch reflections, mixing reflection loss with absorption in a transmission-based estimate. The four-row comparison below was not produced by running that full heterogeneous mapping. It instead tests the `psi=0` pure-cortical-bone limit as a single homogeneous medium, which removes impedance-mismatch reflections from the absorption-calibration check.

For pure cortical-bone-equivalent parameters (c=2800 m/s, rho=1900 kg/m³, alpha_prefactor=21.5, y=1.1, f0=0.5 MHz):

Read the table: each row applies a different combination of dispersion mode (default vs. `no_dispersion`) and whether the coefficient was frequency-matched, holding the medium parameters fixed. Comparing rows shows the individual and combined effect of each fix.

| Mode | alpha_coeff used | Simulated attenuation (dB/cm) | Error vs. literature target | R² |
|---|---:|---:|---:|---:|
| Default dispersive | 1.867490 | 0.775 | −16.96% | ≈1 |
| Default dispersive, corrected coefficient | 2.001526 | 0.823 | −11.89% | ≈1 |
| `no_dispersion`, uncorrected coefficient | 1.867490 | 0.867 | −7.19% | ≈1 |
| `no_dispersion`, corrected coefficient | 2.001526 | 0.929 | −0.54% | ≈1 |

The 0.08686 conversion is correct (converting Np/m to dB/cm: 1 Np ≈ 8.686 dB and 1 m = 100 cm). Omitting frequency matching is a major error source. Correct matching with default dispersion still leaves about −11.89%; `no_dispersion` plus matching reduces it to about −0.54%. The −11.89% differs from the −13.15% reported for y=1.1 in Part 1 because the medium parameters differ (c=3100 m/s, rho=2200 kg/m³ in Part 1 versus c=2800 m/s, rho=1900 kg/m³ here). These values apply only to this synthetic transfer-error example, not to a validated CT-to-acoustic-property model.

<img src="figures/fig6_y11_frequency_scan.png" width="560" alt="Frequency scan for y=1.1 using no-dispersion mode">

Figure 6. Line chart. X-axis: frequency (MHz), 0.2 to 1.0. Y-axis: attenuation (dB/cm). Dashed line: literature y=1 target. Solid line: y=1.1 simulated with no_dispersion mode, using the matched coefficient. The two lines meet at f0=0.5 MHz but diverge at the band edges, showing that no_dispersion fixes the calibration point but not the full-band power-law mismatch.

## Conclusion

Literature attenuation coefficients should be matched to k-Wave's power-law parameterization at the intended operating frequency. This algebraic correction is exact for any exponent k-Wave supports (0 < y ≤ 3), with y=1 handled as a special case. Numerical attenuation near y=1 still requires an independent check. In this narrowband 1-D setup, `no_dispersion` provides a practical alternative. The trade-off is altered phase behavior.

This study does not determine which attenuation exponent is physically appropriate for skull. It examines how a literature coefficient is transferred into k-Wave for a given exponent, and how the resulting peak-pressure attenuation estimate behaves in the simulation.

## Reproduction

### 1. Prerequisites

Install MATLAB R2026a and k-Wave Toolbox 1.4.1. Also install Python 3.11 or later. MATLAB performs the k-Wave simulation; Python only launches MATLAB in batch mode and handles result/figure orchestration. The Python driver does not use the Python k-Wave interface or `matlab.engine`. This project was tested and validated only on MATLAB R2026a. k-Wave 1.4.1 itself may have broader MATLAB version compatibility, but this project has not verified behavior on other MATLAB versions.

### 2. Configure the environment

Make sure `matlab` is available on your system PATH and that the k-Wave Toolbox is discoverable by the MATLAB scripts. From the repository root, install the Python plotting dependencies:

```bash
pip install -r requirements.txt
```

### 3. Run the benchmark

Launch the complete 1-D experiment from the repository root:

```bash
python src/run_experiments.py
```

The command calls `matlab -batch`, writes JSON summaries to `results/`, and regenerates six figures under `figures/`. No GPU is required.

If MATLAB is not available, the committed JSON summaries are sufficient to redraw the figures:

```bash
python src/run_experiments.py --plot-only
```

This requires only `numpy` and `matplotlib` and is the cheapest way to check that the committed results and the published figures are consistent. It cannot regenerate the underlying simulations.

Runtime depends on local MATLAB and k-Wave performance and has not been measured systematically for this project.

**Tested environment and outputs.** The validated setup is MATLAB R2026a, k-Wave Toolbox 1.4.1, and Python 3.11/3.12 on Windows; no GPU is required for the 1-D test. The run writes JSON summaries in `results/` and six figures (`fig1`–`fig6`) under `figures/`.

**Data provenance.** Not every committed JSON is produced by `run_experiments.py`; the gaps are listed here so that no result file is unattributed:

| File | Produced by |
|---|---|
| `results/results.json`, `results/matlab_results.json`, `results/results_n_sweep.json` | `src/run_experiments.py` |
| `results/controls_y098.json` | `src/matlab/run_controls.m` (called by the driver) |
| `results/near_y_frequency_y*.json` | `src/matlab/run_near_y_frequency_scans.m` |
| `results/cuda_y1_backend_control_20260914.json` | `src/matlab/validate_kwave_cuda_y1_*.m` |
| `results/singularity_scan.json`, `results/y11_no_dispersion.json` | Not produced by any committed script. Both were generated by an earlier one-off run that is not included in this repository; they can be re-plotted but not recomputed from the code here. |
| `results/tan_scaling_check.json`, `results/grid_resolution_check.json` | Not produced by any committed script. Exploratory coefficient- and grid-scaling checks for y=1.05; re-plotting is not implemented and the raw values are reported only in these files. |

### 4. Optional CUDA backend control

The CUDA control is separate from the core reproduction and requires a compatible standalone k-Wave CUDA binary. Set `TFUS_KWAVE_ROOT`, `TFUS_KWAVE_BINARY_PATH`, and `TFUS_KWAVE_BINARY_NAME`, then run `src/matlab/validate_kwave_cuda_y1_nearby.m` and `src/matlab/validate_kwave_cuda_y1_timestep.m` from MATLAB. The committed `results/cuda_y1_backend_control_20260914.json` is a path-sanitized summary of the tested run; these controls do not read production data.

## Limitations

- The attenuation-calibration benchmark is a homogeneous 1-D medium. A separate small homogeneous 3-D control tests backend behavior only; neither experiment represents 3-D skull geometry.
- The y≈1 sensitivity boundaries (y=0.80–1.32) depend on this grid,
  PML, and fitting interval. They are not universal constants. The edges are bracketed only by the sampled exponents, so the underlying boundaries lie somewhere in (0.78, 0.80) and (1.32, 1.35).
- In every 1-D run, including the Part 1 results, the source sits inside the internal PML (grid index 8 with a 20-point PML). This is not standard k-Wave practice and is not ruled out as a contributor to the residual few-percent bias seen even at y=2, or to the non-proportional scaling of the attenuation estimate with `alpha_coeff`. Moving the source into the physical interior and lengthening the domain is a required follow-up.
- The attenuation estimate is not proportional to `alpha_coeff` near y=1 (see "Amplitude dependence" in Part 2), so the percentages quoted there are conditional on the coefficient used. The same scaling check has not been run at other exponents.
- `no_dispersion` removes frequency-dependent phase coupling. This
  has little effect for narrowband or CW sources, since there is no
  broadband phase content to distort. It can distort phase for short
  broadband pulses.
- The tan(πy/2) pole is dimension-independent, but the small 3-D control does not establish its numerical severity for other grids, geometries, or binaries.
- Attenuation is estimated from peak pressure rather than a fully
  frequency-resolved spectral attenuation measurement.
- The current experiments do not completely separate dispersion/model
  effects from numerical discretization effects.
- This calibration addresses only the parameter transfer between a literature linear attenuation law and k-Wave's absorption parameters (`alpha_coeff` and `alpha_power`). In real skull, especially trabecular bone, scattering and longitudinal-to-shear mode conversion can also contribute substantially to acoustic energy loss. Those mechanisms are not represented by `alpha_coeff`. Correct absorption calibration therefore does not guarantee agreement with experimentally reported total transmission loss when scattering and mode conversion are present.

**Project structure.** `src/run_experiments.py` orchestrates the core 1-D run, and `src/matlab/run_kwave_experiments.m` performs that k-Wave simulation. The isolated backend controls are `src/matlab/validate_kwave_cuda_y1_nearby.m` and `src/matlab/validate_kwave_cuda_y1_timestep.m`. `results/` stores summaries, and `figures/` stores plots.

## Author

Yue Wang ([Yuewang04](https://github.com/Yuewang04))

## Citation

Machine-readable metadata is provided in [`CITATION.cff`](CITATION.cff). A DOI via Zenodo will be added here once available.

## References

- Aubry, J.-F. et al. “Benchmark problems for transcranial ultrasound simulation: Intercomparison of compressional wave models.” *JASA* 152(2), 1003–1019 (2022). DOI: [10.1121/10.0013426](https://doi.org/10.1121/10.0013426)
- Treeby, B. E. & Cox, B. T. “Modeling power law absorption and dispersion for acoustic propagation using the fractional Laplacian.” *JASA* 127(5), 2741–2748 (2010). DOI: [10.1121/1.3377056](https://doi.org/10.1121/1.3377056)
- Treeby, B. E. & Cox, B. T. “Modeling power law absorption and dispersion in viscoelastic solids using a split-field and the fractional Laplacian.” *JASA* 136(4), 1490–1500 (2014). DOI: [10.1121/1.4894790](https://doi.org/10.1121/1.4894790)
- Singh, A., Jiménez-Gambín, S. & Konofagou, E. E. “An all-ultrasound cranial imaging method to establish the relationship between cranial FUS incidence angle and transcranial attenuation in non-human primates in 3D.” *Scientific Reports* 14, 1488 (2024). DOI: [10.1038/s41598-024-51623-5](https://doi.org/10.1038/s41598-024-51623-5)
- Bamber, J. C. “Attenuation and Absorption.” In C. R. Hill, J. C. Bamber & G. R. ter Haar (Eds.), *Physical Principles of Medical Ultrasonics*, 2nd ed., Chapter 4. John Wiley & Sons (2004). DOI: [10.1002/0470093978.ch4](https://doi.org/10.1002/0470093978.ch4)
- k-Wave documentation, “The k-Wave Manual / power-law absorption and dispersion.” [k-wave.org](https://www.k-wave.org/documentation.php)

## Feedback and license

Questions, corrections, and reproducibility reports are welcome through [GitHub Issues](https://github.com/Yuewang04/kwave-power-law-attenuation-calibration/issues). Licensed under the MIT License; see [LICENSE](LICENSE).

