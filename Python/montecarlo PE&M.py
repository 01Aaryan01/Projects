# scenario      probability     cost_impact
# risk 1        0.15            5
# risk 2        0.2             10
# risk 3        0.1             30
# risk 4        0.25            15
# risk 5        0.15            40
# no risk       0.15            0

# Total project cost with no risk included = 100M
import random
import matplotlib.pyplot as plt
import numpy as np

sim = []

N = int(input("Enter number of simulations: "))

for i in range(N):
    rand = random.random()
    if rand < 0.15:
        cost_impact = 5
    elif rand < 0.35:
        cost_impact = 10
    elif rand < 0.45:
        cost_impact = 30
    elif rand < 0.7:
        cost_impact = 15
    elif rand < 0.85:
        cost_impact = 40
    else:
        cost_impact = 0
    total_cost = 100 + cost_impact
    sim.append(total_cost)

sim.sort()
p90_index = int(0.9 * N) - 1
p90_cost = sim[p90_index]
p80_index = int(0.8 * N) - 1
p80_cost = sim[p80_index]

mean = sum(sim) / N
variance = sum((x - mean) ** 2 for x in sim) / N
std_dev = variance ** 0.5

print(f"P90 Cost: ${p90_cost:,.2f}")
print(f"P80 Cost: ${p80_cost:,.2f}")
print(f"Mean Cost: ${mean:,.2f}")
print(f"Standard Deviation: ${std_dev:,.2f}")
# plt.hist(sim, bins=30, edgecolor='black')
# plt.title('Monte Carlo Simulation of Project Cost')
# plt.xlabel('Total Project Cost ($)')
# plt.ylabel('Frequency')
# plt.axvline(p90_cost, color='r', linestyle='dashed', linewidth=1, label='P90 Cost')
# plt.axvline(p80_cost, color='g', linestyle='dashed', linewidth=1, label='P80 Cost')
# plt.legend()
# plt.show()

x = np.linspace(0, N, N)
y = np.array(sim)
plt.plot(x, y)
plt.title('Monte Carlo Simulation of Project Cost')
plt.xlabel('Simulation Number')
plt.ylabel('Total Project Cost ($)')
plt.axhline(p90_cost, color='r', linestyle='dashed', linewidth=1, label='P90 Cost')
plt.axhline(p80_cost, color='g', linestyle='dashed', linewidth=1, label='P80 Cost')
plt.axhline(mean, color='b', linestyle='dashed', linewidth=1, label='Mean Cost')
# plt.axhline(mean + std_dev, color='y', linestyle='dashed', linewidth=1, label='Mean + 1 Std Dev')
# plt.axhline(mean - std_dev, color='y', linestyle='dashed', linewidth=1, label='Mean - 1 Std Dev')
plt.legend()
plt.grid()
plt.show()



