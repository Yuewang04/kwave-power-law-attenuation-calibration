# Near-y=1 Experiment Protocol & Completed Controls

## Objective

Diagnose why k-Wave's default dispersive mode produces large attenuation errors and field instabilities near $y=1$, and reconcile 1-D and 3-D behavior using exact fractional-Laplacian PDE theory.

## Completed Standard Setup (v1.1.0)

- **Solver:** MATLAB R2026a, k-Wave 1.4.1 (`kspaceFirstOrder1D`, `'DataCast', 'double'`).
- **Medium:** `c = 3100 m/s`, `rho = 2200 kg/m^3`, domain `256 mm` (`Nx = 512`, `dx = 0.5 mm`), `PML = 20` grid points.
- **Source & Duration:** 18-cycle toneburst placed strictly in the physical interior at grid index **35** (15 cells inside the physical domain, outside the 20-cell PML), simulated for **180 µs** so even dispersively slowed waves (`cp ≈ 1767 m/s` at `y = 0.98`) completely clear the farthest sensor.
- **Sensors:** 6 attenuation probes at `26, 46, 66, 86, 106, 126 mm` from the source, plus an adjacent `+1 dx` (`0.5 mm`) probe for unaliased phase-speed extraction. Recorded fields: `{'p', 'p_max', 'p_min'}` (`max(|p_max|, |p_min|)`).
- **Frequencies:** `0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0 MHz`.

## Completed Control Checks

1. **Timestep refinement (`dt`, `dt/2`, `dt/4`) at `y = 0.98` and `y = 0.99`** (`src/matlab/run_controls.m` -> `results/controls_y098.json`):
   - `y = 0.98` remains stable (`+190.45%`, `+190.18%`, `+190.15%`, `R² = 1.0000`) and matches spectral FFT estimates (`+171–173%`).
   - `y = 0.99` unconditionally diverges across all timesteps (`max|p| ~ 10^301 Pa`, negative-going elliptic blow-up because effective stiffness $M(k_0) = -0.3505 < 0$).
2. **Grid refinement (`dx = 0.5, 0.25, 0.125 mm`) at `y = 1.05`** (`results/grid_resolution_check.json`):
   - Errors remain `-22.44%`, `-22.36%`, `-22.27%`, matching the Treeby & Cox (2014) analytical prediction (`-21.65%`).
3. **Coefficient scaling (`0.5x, 1x, 2x`) at `y = 1.05` vs `y = 2.0`** (`results/tan_scaling_check.json`):
   - Nonlinear scaling at `y = 1.05` (`-12.50%, -22.44%, -36.65%`) tracks the analytical dispersion denominator `1 - (y+1)*a0*c0*tan(pi*y/2)*w^(y-1)` while spatial decay stays strictly exponential (`R² >= 0.999988`).
4. **Full-band `y = 1.0` with `no_dispersion`** (`results/y1_no_dispersion_fullband.json`):
   - Error stays within `-0.16%` to `-1.59%` across `0.2–1.0 MHz` with unshifted `cp = 3100–3105 m/s`.
5. **Audit of k-Wave built-in `fitPowerLawParams`** (`results/fit_power_law_params_audit.json`):
   - Recovers target attenuation amplitude near `y = 1` (`-0.78%` to `-1.41%` for `y_desired in [1.00, 1.10]`), but introduces `+6.5%` to `+103.7%` superluminal phase-speed distortion (`cp = 3303–6316 m/s`), and still diverges for `y_desired = 0.99`.

Reproduce the entire 1-D suite in MATLAB via `run_all_benchmarks` (or `python src/run_experiments.py`), or redraw all figures without MATLAB via `python src/run_experiments.py --plot-only` and `python src/plot_near_y1_frequency_scan.py`.
