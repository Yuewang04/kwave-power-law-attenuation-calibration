# Near-y=1 frequency-scan experiment

## Objective

Test whether the large attenuation-estimation error observed at 500 kHz is a frequency-local effect or a systematic near-y=1 behavior.

## Completed matrix

- Solver: MATLAB k-Wave 1-D, `kspaceFirstOrder1D`.
- Medium: `c=3100 m/s`, `rho=2200 kg/m^3`, domain `256 mm`, baseline `dx=0.5 mm`, `PML=20` grid points.
- Source: 18-cycle Gaussian toneburst, source at grid index 8.
- Sensors: 26, 46, 66, 86, 106, and 126 mm from the source.
- Frequencies: 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, and 1.0 MHz.
- Exponents added in this run: `y=0.98, 0.99, 1.01, 1.02`.
- The coefficient is rescaled at `f0=0.5 MHz` for every exponent. The reported observable is the linear-fit attenuation from `ln(p_max)` versus distance.

Existing saved curves for `y=0.9, 1.0, 1.05, 1.1` use the same frequency grid and configuration. Exact `y=1` uses `alpha_mode='no_dispersion'`; the near-y=1 cases use the default dispersive mode.

## Required controls before interpretation

1. Repeat `y=0.98, 0.99, 1.01, 1.02` at `CFL=0.05` and `CFL=0.025`; the baseline currently emits a time-step warning for some near-y=1 cases.
2. Repeat at `dx=0.25 mm` with the time step reduced consistently.
3. Save raw sensor time series and compare peak-pressure, gated FFT, and cross-spectral attenuation estimates.
4. Add matched `alpha_coeff=0` runs to remove source/grid/PML transfer effects.
5. Keep `y=1` as a MATLAB `no_dispersion` reference and test CUDA behavior separately; do not infer GPU portability from the CPU result.

The source lies inside the internal PML. At 0.2 MHz the 18-cycle burst lasts approximately 90 microseconds, matching the requested record duration, so distant sensors cannot record its complete tail. Move the source into the physical interior and lengthen the record in controlled reruns.

The current matrix is suitable for locating frequency dependence. It is not by itself a proof of intrinsic material attenuation or a universal safe distance from `y=1`. Errors use the linear literature target; away from the fixed matching frequency of 0.5 MHz they include the expected power-law mismatch `(f/f0)^(y-1)-1`. Values below -100% correspond to negative fitted attenuation, not negative material absorption.

Reproduce with `run_near_y_frequency_scans` after adding `src/matlab` to the MATLAB path, then run `python src/plot_near_y1_frequency_scan.py`.
