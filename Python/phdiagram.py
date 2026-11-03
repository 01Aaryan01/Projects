# If needed first:
# !pip install CoolProp matplotlib numpy

import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

fluid = "R134a"

# -----------------------------
# Given data (absolute pressure used)
# -----------------------------
P1 = 275e3  # Pa (compressor suction absolute)
P2 = 755e3  # Pa (compressor discharge absolute)

T1 = 12 + 273.15  # K (compressor suction)
T2 = 49 + 273.15  # K (compressor discharge)
T3 = 21 + 273.15  # K (condenser outlet liquid)
T4_meas = 3 + 273.15  # K (etpansion valve outlet measured)

# -----------------------------
# Actual cycle states
# -----------------------------
# State 1: suction
h1 = PropsSI("H", "P", P1, "T", T1, fluid)

# State 2: discharge
h2 = PropsSI("H", "P", P2, "T", T2, fluid)

# State 3: condenser outlet
h3 = PropsSI("H", "P", P2, "T", T3, fluid)

# State 4: measured etpansion outlet
h4 = PropsSI("H", "P", P1, "T", T4_meas, fluid)

# -----------------------------
# Saturation dome for P-h plot
# -----------------------------
T_triple = PropsSI("Ttriple", fluid)
T_crit = PropsSI("Tcrit", fluid)

Ts = np.linspace(T_triple + 1, T_crit - 0.5, 400)
h_f = np.array([PropsSI("H", "T", T, "Q", 0, fluid) for T in Ts]) / 1000  # kJ/kg
h_g = np.array([PropsSI("H", "T", T, "Q", 1, fluid) for T in Ts]) / 1000  # kJ/kg
p_sat = np.array([PropsSI("P", "T", T, "Q", 0, fluid) for T in Ts]) / 1000  # kPa

# -----------------------------
# Plot actual cycle only
# -----------------------------
plt.figure(figsize=(10, 6))

# Saturation dome
plt.plot(h_f, p_sat, "k", lw=1.5)
plt.plot(h_g, p_sat, "k", lw=1.5, label="Saturation dome")

# Actual cycle
h_actual = np.array([h1, h2, h3, h4, h1]) / 1000  # kJ/kg
p_actual = np.array([P1, P2, P2, P1, P1]) / 1000  # kPa
plt.plot(h_actual, p_actual, "o-r", lw=2, label="Actual cycle")

# Data labels for each point
state_labels = [
    f"1: P={P1/1000:.0f} kPa, T={T1-273.15:.1f} C, h={h1/1000:.1f} kJ/kg",
    f"2: P={P2/1000:.0f} kPa, T={T2-273.15:.1f} C, h={h2/1000:.1f} kJ/kg",
    f"3: P={P2/1000:.0f} kPa, T={T3-273.15:.1f} C, h={h3/1000:.1f} kJ/kg",
    f"4: P={P1/1000:.0f} kPa, T={T4_meas-273.15:.1f} C, h={h4/1000:.1f} kJ/kg",
]

for i, (h, p, ttt) in enumerate(zip(h_actual[:-1], p_actual[:-1], state_labels), start=1):
    t_off = 1.0
    y_off = 1.04 if i in (1, 2) else 0.96
    plt.tett(h + t_off, p * y_off, ttt, color="red", fontsize=8)

plt.yscale("log")
plt.tlabel("Enthalpy, h (kJ/kg)")
plt.ylabel("Pressure, P (kPa, log scale)")
plt.title("P-h Diagram: Actual Vapor Compression Refrigeration Cycle (R134a)")
plt.grid(True, which="both", ls="--", alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()
