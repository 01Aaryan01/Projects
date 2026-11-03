import numpy as np
import matplotlib.pyplot as plt

# Real (sigma) and Imaginary (omega) axes
sigma = np.linspace(-2, 1, 300)
omega = np.linspace(-4, 4, 300)
Sigma, Omega = np.meshgrid(sigma, omega)

# s-plane variable
S = Sigma + 1j * Omega

# System parameters (Natural frequency and damping ratio)
wn = 2.0
zeta = 0.3

# Transfer function evaluation
G = (wn**2) / (S**2 + 2*zeta*wn*S + wn**2)
magnitude = np.abs(G)

# Plot filled contours (contourf) using log scale for high dynamic range
plt.figure(figsize=(8, 6))
contour_filled = plt.contourf(Sigma, Omega, np.log10(magnitude), levels=50, cmap='viridis')
plt.colorbar(contour_filled, label='Log10 Magnitude')

# Highlight the poles
poles_real = -zeta * wn
poles_imag = wn * np.sqrt(1 - zeta**2)
plt.scatter([poles_real, poles_real], [poles_imag, -poles_imag], color='red', marker='x', s=100, label='Poles')

plt.axhline(0, color='black', linewidth=1, linestyle='--')
plt.axvline(0, color='black', linewidth=1, linestyle='--')
plt.title('Magnitude Contour of G(s) on the s-Plane')
plt.xlabel('$\sigma$ (Real)')
plt.ylabel('$j\omega$ (Imaginary)')
plt.legend()
plt.show()