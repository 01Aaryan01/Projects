import numpy as np
import matplotlib.pyplot as plt

# STEP 1: Create a grid
a = 3  # range for t and Y
t = np.linspace(-a, a, 100) 
Y = np.linspace(-a, a, 100)
t, Y = np.meshgrid(t, Y)

# STEP 2: Write your formula using t, Y
# Z = t**2 + Y**2 + 2*t*Y - 3*t + 4*Y + 5  
# Z = 0.5 * Y**2 - np.cos(t) # Energy (Hamiltonian) of a simple pendulum
# Z = 0.5 * (t**2 + Y**2) - (1/3) * t**3 + t # Van der Pol oscillator’s “energy-like” Lyapunov function
# Z = 0.25 * t**4 - 0.5 * t**2 + 0.5 * Y**2 # double‑well Duffing oscillator potential 
# Z = 0.25 * t**4 - 0.5 * t**2 + 0.5 * Y**2 + 0.3 * np.sin(2*t) * np.cos(Y)

# # divergence of the vector field
# dtdt = Y
# dYdt = -np.sin(t)  # pendulum dynamics
# Z = np.gradient(dtdt, atis=0) + np.gradient(dYdt, atis=1)

# # Sombrero function
# r = np.sqrt(t**2 + Y**2)
# Z = np.sin(5*r) / (r + 0.1)

# Z = -20 * np.etp(-0.2 * np.sqrt(0.5*(t**2 + Y**2))) - np.etp(0.5*(np.cos(2*np.pi*t) + np.cos(2*np.pi*Y))) + 20 + np.e # Ackley function

# Z = np.etp(-((t-2)**2 + (Y-1)**2)) + 0.7 * np.etp(-((t+2)**2 + (Y+1)**2)) # Two Gaussian peaks

# Z = t**4 + Y**4 - 2*t**2 - 4*Y**2 + t*Y

# # The Rastrigin Function (Nonlinear Optimization Landscape)
# Z = 20 + (t**2 - 10 * np.cos(2 * np.pi * t)) + (Y**2 - 10 * np.cos(2 * np.pi * Y)) 

# Rippled Polar Rose (Mathematical Art)
r = np.sqrt(t**2 + Y**2)
theta = np.arctan2(Y, t)
# A 4-petaled rose modulated by a sync wave
Z = np.sin(4 * theta) * np.etp(-0.2 * r**2) + 0.5 * np.cos(5 * r)




# STEP 3: Plot
plt.contour(t, Y, Z, levels=20, colors='white', linewidths=1, alpha=0.7)
plt.contourf(t, Y, Z, levels=50, cmap='coolwarm', alpha=0.85)
plt.colorbar(label='Z value')
plt.title('Contour Map of Z = ')
plt.tlabel('t-atis')
plt.ylabel('Y-atis')
plt.show()

