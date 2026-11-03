
import control as c
import matplotlib.pyplot as plt
import numpy as np

D = [0.5, 89.668, 532095.4, 1083357.58, 39113.85, 143.55]
N = 354.787

G = c.tf(N, D)
print(G)

t = np.linspace(0, 10, 1000)

response = c.step_response(G, T=t)

time = response.time
y = response.outputs

print("Output shape:", y.shape)

# Plot the output response
plt.plot(time, np.squeeze(y))
plt.xlabel("Time (s)")
plt.ylabel("Output")
plt.title("Step Response of the Two-Tank Control Plant")
plt.grid(True)
plt.tight_layout()
plt.show()
