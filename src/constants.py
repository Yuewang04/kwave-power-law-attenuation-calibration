"""Physical constants and unit conversions for the focused experiment."""

c_water = 1500.0
c_bone = 3100.0
rho_water = 1000.0
rho_bone = 2200.0
alpha_min_bone_Np = 21.5  # Np/MHz/m, literature y=1 coefficient
Np2dB = 8.686
alpha_lit_dB_cm_MHz = alpha_min_bone_Np * Np2dB / 100.0
f0_MHz = 0.5
alpha_power = 2
alpha_coeff_uncorrected = alpha_lit_dB_cm_MHz
alpha_coeff_corrected = alpha_lit_dB_cm_MHz / f0_MHz
target_dB_cm_at_f0 = alpha_lit_dB_cm_MHz * f0_MHz
