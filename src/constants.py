"""Physical constants, exact unit conversions, and analytical fractional-Laplacian predictions."""
import math

c_water = 1500.0
c_bone = 3100.0
rho_water = 1000.0
rho_bone = 2200.0
alpha_min_bone_Np = 21.5  # Np / (MHz * m), literature y=1 coefficient
Np2dB = 20.0 * math.log10(math.e)  # Exact ~8.685889638, matching k-Wave's db2neper.m
alpha_lit_dB_cm_MHz = alpha_min_bone_Np * Np2dB / 100.0
f0_MHz = 0.5
alpha_power = 2
alpha_coeff_uncorrected = alpha_lit_dB_cm_MHz
alpha_coeff_corrected = alpha_lit_dB_cm_MHz / f0_MHz
target_dB_cm_at_f0 = alpha_lit_dB_cm_MHz * f0_MHz


def db2neper(alpha_db_cm_mhzy: float, y: float) -> float:
    """Convert dB/(cm*MHz^y) to Np/((rad/s)^y * m), matching k-Wave db2neper.m."""
    return 100.0 * float(alpha_db_cm_mhzy) * ((1e-6 / (2.0 * math.pi)) ** float(y)) / Np2dB


def neper2db(alpha_np: float, y: float) -> float:
    """Convert Np/((rad/s)^y * m) back to dB/(cm*MHz^y), matching k-Wave neper2db.m."""
    return float(alpha_np) * Np2dB / (100.0 * ((1e-6 / (2.0 * math.pi)) ** float(y)))


def matched_alpha_coeff(alpha_lit: float, y: float, f0: float = f0_MHz) -> float:
    """Match a linear (y=1) literature coefficient at reference frequency f0 (MHz)."""
    return float(alpha_lit) * (float(f0) ** (1.0 - float(y)))


def treeby_stiffness_multiplier(f_MHz: float, alpha_coeff: float, y: float, c0: float) -> float:
    """Return effective stiffness multiplier M(k0) = 1 - eta * k0^(y-1) in k-Wave's state equation.

    When M <= 0 (e.g. y < 1 sufficiently close to 1), c_eff^2 <= 0 and the
    hyperbolic wave equation becomes Hadamard ill-posed (elliptic blow-up).
    """
    if abs(y - 1.0) < 1e-12:
        return float("nan")
    a0_np = db2neper(alpha_coeff, y)
    w = 2.0 * math.pi * float(f_MHz) * 1e6
    k0 = w / float(c0)
    eta = 2.0 * a0_np * (float(c0) ** float(y)) * math.tan(math.pi * float(y) / 2.0)
    return 1.0 - eta * (k0 ** (float(y) - 1.0))


def treeby_actual_attenuation_db_cm(
    f_MHz: float, alpha_coeff: float, y: float, c0: float, no_dispersion: bool = False
) -> float:
    """Return exact analytical attenuation (dB/cm) of k-Wave's first-order fractional PDE.

    From Treeby & Cox (2014) / k-Wave fitPowerLawParams.m:
      alpha_actual(w) = a0*w^y / (1 - (y+1)*a0*c0*tan(pi*y/2)*w^(y-1))
    Returns NaN if the dispersion pole drives M(k0) <= 0 or denominator <= 0.
    """
    if no_dispersion:
        return float(alpha_coeff) * (float(f_MHz) ** float(y))
    if abs(y - 1.0) < 1e-12:
        return float("nan")
    a0_np = db2neper(alpha_coeff, y)
    w = 2.0 * math.pi * float(f_MHz) * 1e6
    tan_val = math.tan(math.pi * float(y) / 2.0)
    denom = 1.0 - (float(y) + 1.0) * a0_np * float(c0) * tan_val * (w ** (float(y) - 1.0))
    m_k0 = treeby_stiffness_multiplier(f_MHz, alpha_coeff, y, c0)
    if denom <= 0.0 or m_k0 <= 0.0:
        return float("nan")
    alpha_np_m = (a0_np * (w ** float(y))) / denom
    return alpha_np_m * 0.01 * Np2dB


def treeby_phase_speed_m_s(
    f_MHz: float, alpha_coeff: float, y: float, c0: float, no_dispersion: bool = False
) -> float:
    """Return low-loss analytical phase speed c_eff = c0 * sqrt(M(k0)) in k-Wave."""
    if no_dispersion:
        return float(c0)
    m_k0 = treeby_stiffness_multiplier(f_MHz, alpha_coeff, y, c0)
    if not math.isfinite(m_k0) or m_k0 <= 0.0:
        return float("nan")
    return float(c0) * math.sqrt(m_k0)
