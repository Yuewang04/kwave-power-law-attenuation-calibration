"""Small, deterministic 1-D attenuation estimator fixture and analytical verification tests.

!! WARNING -- NOT A K-WAVE VALIDATION TOOL !!

This module is an *analytical* surrogate, not a k-Wave simulation. It builds
``exp(-alpha * x)`` directly, so it returns the target attenuation by
construction: every "calibration error" computed from it is zero to machine
precision. It must never be used to check whether k-Wave reproduces a
converted coefficient, and it is NOT imported by ``run_experiments.py``.

It is retained as a unit-test fixture for the log-fit estimator and the
exact fractional-Laplacian PDE formulas in ``constants.py``. The real 1-D
k-Wave experiments live in ``src/matlab/``.
"""
import math
import numpy as np

from constants import (
    Np2dB,
    alpha_lit_dB_cm_MHz,
    c_bone,
    db2neper,
    matched_alpha_coeff,
    neper2db,
    treeby_actual_attenuation_db_cm,
    treeby_stiffness_multiplier,
)


def simulate_attenuation(f_MHz, alpha_coeff, alpha_power=2,
                         distances_m=None, cycles=18):
    """Return fitted attenuation (dB/cm) and diagnostic arrays.

    WARNING: this is a closed-form plane-wave decay, not a simulation.
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


def run_self_tests():
    """Verify log-fit recovery, unit conversions, and analytical PDE regime signs."""
    # 1. Closed-form log-fit recovery
    a0_y2 = matched_alpha_coeff(alpha_lit_dB_cm_MHz, 2.0, 0.5)
    out = simulate_attenuation(0.5, a0_y2, alpha_power=2.0)
    target = alpha_lit_dB_cm_MHz * 0.5
    assert abs(out["attenuation_dB_cm"] - target) < 1e-12, "Log-fit recovery failed"

    # 2. db2neper <-> neper2db round-trip
    for y in (0.5, 0.99, 1.0, 1.05, 2.0):
        rt = neper2db(db2neper(alpha_lit_dB_cm_MHz, y), y)
        assert abs(rt - alpha_lit_dB_cm_MHz) < 1e-12, f"Unit roundtrip failed at y={y}"

    # 3. Treeby analytical formula checks
    att_y2 = treeby_actual_attenuation_db_cm(0.5, a0_y2, 2.0, c_bone)
    assert abs(att_y2 - target) < 1e-12, "y=2 tan(pi)=0 must give exact target"

    a0_101 = matched_alpha_coeff(alpha_lit_dB_cm_MHz, 1.01, 0.5)
    att_101 = treeby_actual_attenuation_db_cm(0.5, a0_101, 1.01, c_bone)
    err_101 = (att_101 - target) / target * 100.0
    assert -60.0 < err_101 < -55.0, f"Expected ~-57.6% at y=1.01, got {err_101:.2f}%"

    # 4. Hadamard instability check: M(k0) < 0 at y=0.99
    a0_099 = matched_alpha_coeff(alpha_lit_dB_cm_MHz, 0.99, 0.5)
    m_099 = treeby_stiffness_multiplier(0.5, a0_099, 0.99, c_bone)
    assert m_099 < 0.0, f"Expected negative effective stiffness at y=0.99, got {m_099}"
    assert math.isnan(treeby_actual_attenuation_db_cm(0.5, a0_099, 0.99, c_bone))
    print("[✓] All analytical estimator & PDE regime self-tests passed.")


if __name__ == "__main__":
    run_self_tests()
