import numpy as np
import matplotlib.pyplot as plt

# Image resolution & bounds
width, height = 800, 800
xmin, xmax = -2, 2
ymin, ymax = -2, 2
max_iter = 50
tolerance = 1e-6

# Create grid of complex numbers
x = np.linspace(xmin, xmax, width)
y = np.linspace(ymin, ymax, height)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y

# Define the function and derivative
def f(z):
    return z**3 - 1

def f_prime(z):
    return 3*z**2

# Known roots
roots = np.array([
    1 + 0j,
    -0.5 + np.sqrt(3)/2 * 1j,
    -0.5 - np.sqrt(3)/2 * 1j
])

# Store root index & iterations
root_index = np.zeros(Z.shape, dtype=int)
iterations = np.zeros(Z.shape, dtype=int)

# Newton iteration
for k in range(max_iter):
    Z -= f(Z) / f_prime(Z)
    for i, r in enumerate(roots):
        mask = np.abs(Z - r) < tolerance
        root_index[mask] = i
        iterations[mask] = k

# Normalize colors
colors = (root_index + iterations / max_iter) / len(roots)

# Plot
plt.imshow(colors, extent=[xmin, xmax, ymin, ymax], cmap='twilight')
plt.title("Newton's Fractal for $z^3 - 1 = 0$")
plt.xlabel("Re(z)")
plt.ylabel("Im(z)")
plt.show()
