"""Small, deterministic 1-D attenuation experiment.

!! WARNING -- NOT A VALIDATION TOOL !!

This module is an *analytical* surrogate, not a k-Wave simulation. It builds
``exp(-alpha * x)`` directly, so it returns the target attenuation by
construction: every "calibration error" computed from it is zero to machine
precision. It must never be used to check whether k-Wave reproduces a
converted coefficient, and it is NOT imported by ``run_experiments.py``.

It is retained only as a unit-test fixture for the log-fit estimator. The real
1-D experiment lives in ``src/matlab/run_kwave_experiments.m``.
"""
import numpy as np

from constants import Np2dB


def simulate_attenuation(f_MHz, alpha_coeff, alpha_power=2,
                        distances_m=None, cycles=18):
    """Return fitted attenuation (dB/cm) and diagnostic arrays.

    WARNING: this is a closed-form plane-wave decay, not a simulation. See the
    module docstring.

    A toneburst's peak amplitude decays as exp(-alpha_np_m*x), exactly the
    homogeneous-medium solution represented by k-Wave's power-law model.
    ``cycles`` documents the source setup (18-cycle toneburst); peak extraction
    is represented by the envelope because the medium is uniform and linear.
    """
    if distances_m is None:
        distances_m = np.linspace(0.02, 0.14, 13)
    distances_m = np.asarray(distances_m, dtype=float)
    dB_cm = float(alpha_coeff) * float(f_MHz) ** float(alpha_power)
    alpha_np_m = dB_cm * 100.0 / Np2dB
    amplitudes = np.exp(-alpha_np_m * distances_m)
    slope, intercept = np.polyfit(distances_m, np.log(amplitudes), 1)
    fitted_np_m = -slope
    fitted_dB_cm = fitted_np_m * Np2dB / 100.0
    return {
        "frequency_MHz": float(f_MHz),
        "alpha_coeff": float(alpha_coeff),
        "alpha_power": float(alpha_power),
        "cycles": int(cycles),
        "attenuation_dB_cm": float(fitted_dB_cm),
        "distances_m": distances_m.tolist(),
        "p_max": amplitudes.tolist(),
        "fit_intercept": float(intercept),
        "fit_slope_Np_per_m": float(slope),
    }
