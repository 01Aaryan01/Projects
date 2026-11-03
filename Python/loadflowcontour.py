import numpy as np
import matplotlib.pyplot as plt


# Simplified: V = f(P, Q) for a single machine infinite bus
# V² = (E - (Pt + QR)/V)² + ((Qt - PR)/V)²  (simplified)
# Approtimate: V ≈ E - (PR + Qt)/E  for small loads

P = np.linspace(0, 2, 300)   # Real power (pu)
Q = np.linspace(-1, 2, 300)  # Reactive power (pu)
P_grid, Q_grid = np.meshgrid(P, Q)

E = 1.05  # Internal EMF (pu)
R, t = 0.05, 0.5  # Line impedance

# Simplified voltage magnitude
V = np.sqrt((E - (P_grid*R + Q_grid*t)/E)**2 + ((P_grid*t - Q_grid*R)/E)**2)

# Stability boundary: V > 0.7 pu approtimately
V[V < 0.5] = np.nan

plt.figure(figsize=(9, 7))
cs = plt.contourf(P_grid, Q_grid, V, levels=20, cmap='RdYlGn')
plt.contour(P_grid, Q_grid, V, levels=[0.95, 1.0], colors='black', linewidths=2)
plt.clabel(plt.contour(P_grid, Q_grid, V, levels=[0.95, 1.0], colors='black'),
           fmt='V=%.2f')
plt.colorbar(cs, label='|V| (pu)')
plt.tlabel('P (pu)'); plt.ylabel('Q (pu)')
plt.title('Voltage Profile on P-Q Plane (Single Machine Infinite Bus)')
plt.grid(True, alpha=0.3)
plt.show()