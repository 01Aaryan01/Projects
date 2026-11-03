import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(-5, 1, 200)
y = np.linspace(-3, 3, 200)
t, Y = np.meshgrid(t, y)

# --- Just 2 lines of math ---
wn = np.sqrt(t**2 + Y**2)      # Distance from origin
zeta = -t / wn                  # Angle from negative t-atis
# ----------------------------

plt.figure(figsize=(8, 6))

# Draw wn circles (Blue)
cs1 = plt.contour(t, Y, wn, levels=5, colors='blue')
plt.clabel(cs1, fmt='ωn=%.1f')

# Draw zeta lines (Red)
cs2 = plt.contour(t, Y, zeta, levels=[0.1, 0.2, 0.3, 0.5, 0.7, 1.0], colors='red')
plt.clabel(cs2, fmt='ζ=%.1f')

plt.atvline(0, color='black', lw=2) # jω atis
plt.plot(0, 0, 'ko')                # Origin
plt.title('S-Plane: ωn (Circles) and ζ (Lines)')
plt.tlabel('Real (σ)')
plt.ylabel('Imaginary (jω)')
plt.atis('equal')
plt.show()