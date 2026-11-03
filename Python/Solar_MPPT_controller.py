# Solar MPPT Controller Simulation using Single Diode Model

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# # Single diode model of a PV cell
# class PV_cell:
#     def __init__(self):
#         # Physical Constants
#         self.q = 1.602176634e-19      # Charge of an electron (C)
#         self.k = 1.380649e-23         # Boltzmann constant (J/K)
        
#         # Reference Conditions (STC)
#         self.T_c_ref = 25 + 273.15    # Reference temp (K)
#         self.G_ref = 1000             # Reference irradiance (W/m^2)
        
#         # Cell Parameters (From Datasheet or Estimation)
#         self.E_g = 1.1 * self.q       # Bandgap energy (J) for Silicon
#         self.I_sc_ref = 5.0           # Short-circuit current at STC (A)
#         self.K_i = 0.005              # Temp coefficient of Isc (A/K)
#         self.V_oc_ref = 0.6           # Open-circuit voltage at STC (V)
#         self.n = 1.3                  # Diode ideality factor
        
#         # Derived Parameters (Calculated once)
#         self.R_s = 0.01               # Series resistance (Ohms)
#         self.R_sh = 1000              # Shunt resistance (Ohms)
#         self.I_0_ref = (self.I_sc_ref - self.V_oc_ref/self.R_sh) / \
#             (np.exp((self.q * self.V_oc_ref) / (self.n * self.k * self.T_c_ref)) - 1)
    
#     def f1(self, I, V, G, T_c): 
#         V_t = self.n * self.k * T_c / self.q # Thermal voltage
#         I_0 = self.I_0_ref * (T_c/self.T_c_ref)**3 * \
#             np.exp((self.q * self.E_g) / (self.n * self.k) * (1/self.T_c_ref - 1/T_c))
#         I_ph = (self.I_sc_ref + self.K_i * (T_c - self.T_c_ref)) * (G / self.G_ref)
#         I_d = I_0 * (np.exp((V + I * self.R_s) / V_t) - 1)
#         I_sh = (V + I * self.R_s) / self.R_sh 

#         return I - I_ph + I_d + I_sh
    
#     def df1_dI(self, I, V, T_c):
#         V_t = self.n * self.k * T_c / self.q # Thermal voltage
#         I_0 = self.I_0_ref * (T_c/self.T_c_ref)**3 * \
#             np.exp((self.q * self.E_g) / (self.n * self.k) * (1/self.T_c_ref - 1/T_c))
        
#         return 1 + (I_0 * np.exp((V + I * self.R_s) / V_t) * (self.R_s / V_t)) + (self.R_s / self.R_sh)


#     # Newton-Raphson method to solve for current I for given V, G and T_c
#     def solve_I(self, V, G, T_c):
#         I_ph = (self.I_sc_ref + self.K_i * (T_c - self.T_c_ref)) * (G / self.G_ref)
#         I_guess = I_ph  # Initial guess (Maximum possible current is Iph)
#         tolerance = 1e-6
#         max_iter = 100

#         for i in range(max_iter):
#             f_val = self.f1(I_guess, V, G, T_c)
#             df_dI_val = self.df1_dI(I_guess, V, T_c)

#             if df_dI_val == 0:  # Avoid division by zero
#                 break

#             I_new = I_guess - f_val / df_dI_val

#             if abs(I_new - I_guess) < tolerance:
#                 return max(0, I_new)  # Current cannot be negative

#             I_guess = I_new
        
#         return max(0, I_guess)  # Return the last guess if max iterations reached   
    
#     def plot_IV(self, G, T_c):
#         voltages = np.linspace(0, self.V_oc_ref, 100)
#         currents = [self.solve_I(V, G, T_c) for V in voltages]
#         currents_low_G = [self.solve_I(V, G - 500, T_c) for V in voltages]
#         currents_high_G = [self.solve_I(V, G + 500, T_c) for V in voltages]
#         currents_low_T = [self.solve_I(V, G, T_c - 20) for V in voltages]
#         currents_high_T = [self.solve_I(V, G, T_c + 20) for V in voltages]
#         plt.plot(voltages, currents, label='At Standard Conditions')
#         plt.plot(voltages, currents_low_G, label='Lower Irradiance (G-500 W/m^2)')
#         plt.plot(voltages, currents_high_G, label='Higher Irradiance (G+500 W/m^2)')
#         plt.plot(voltages, currents_low_T, label='Lower Temperature (T-20 °C)')
#         plt.plot(voltages, currents_high_T, label='Higher Temperature (T+20 °C)')
#         plt.legend()
#         plt.xlabel('Voltage (V)')
#         plt.ylabel('Current (A)')
#         plt.title(f'I-V Curve at G={G} W/m^2 and T_c={T_c-273.15}°C')
#         plt.grid()
#         plt.show()
    
#     def plot_PV(self, G, T_c):
#         voltages = np.linspace(0, self.V_oc_ref, 100)
#         currents = [self.solve_I(V, G, T_c) for V in voltages]
#         currents_low_G = [self.solve_I(V, G - 500, T_c) for V in voltages]
#         currents_high_G = [self.solve_I(V, G + 500, T_c) for V in voltages]
#         currents_low_T = [self.solve_I(V, G, T_c - 20) for V in voltages]
#         currents_high_T = [self.solve_I(V, G, T_c + 20) for V in voltages]
#         powers = voltages * currents
#         powers_low_G = voltages * currents_low_G
#         powers_high_G = voltages * currents_high_G
#         powers_low_T = voltages * currents_low_T
#         powers_high_T = voltages * currents_high_T
#         plt.plot(voltages, powers, label='At Standard Conditions')
#         plt.plot(voltages, powers_low_G, label='Lower Irradiance (G-500 W/m^2)')
#         plt.plot(voltages, powers_high_G, label='Higher Irradiance (G+500 W/m^2)')
#         plt.plot(voltages, powers_low_T, label='Lower Temperature (T-20 °C)')
#         plt.plot(voltages, powers_high_T, label='Higher Temperature (T+20 °C)')
#         plt.legend()
#         plt.xlabel('Voltage (V)')
#         plt.ylabel('Power (W)')
#         plt.title(f'P-V Curve at G={G} W/m^2 and T_c={T_c-273.15}°C')
#         plt.grid()
#         plt.show()
    
# cell = PV_cell()
# # Plot I-V and P-V curves at standard conditions
# cell.plot_IV(cell.G_ref, cell.T_c_ref)
# cell.plot_PV(cell.G_ref, cell.T_c_ref)


# # PV_array model consisting of multiple cells in series and parallel
# class PV_array: 
#     def __init__(self, N_s, N_p):
#         self.N_s = N_s  # Number of cells in series
#         self.N_p = N_p  # Number of cells in parallel

#         # Physical Constants
#         self.q = 1.602176634e-19      # Charge of an electron (C)
#         self.k = 1.380649e-23         # Boltzmann constant (J/K)
        
#         # Reference Conditions (STC)
#         self.T_c_ref = 25 + 273.15    # Reference temp (K)
#         self.G_ref = 1000             # Reference irradiance (W/m^2)
        
#         # Cell Parameters (From Datasheet or Estimation)
#         self.E_g = 1.1 * self.q       # Bandgap energy (J) for Silicon
#         self.I_sc_ref = 5.0           # Short-circuit current at STC (A)
#         self.K_i = 0.005              # Temp coefficient of Isc (A/K)
#         self.V_oc_ref = 0.6           # Open-circuit voltage at STC (V)
#         self.n = 1.3                  # Diode ideality factor
        
#         # Derived Parameters (Calculated once)
#         self.R_s = 0.01               # Series resistance (Ohms)
#         self.R_sh = 1000              # Shunt resistance (Ohms)
#         self.I_0_ref = (self.I_sc_ref - self.V_oc_ref/self.R_sh) / \
#             (np.exp((self.q * self.V_oc_ref) / (self.n * self.k * self.T_c_ref)) - 1)
    
#     def f1(self, I, V, G, T_c): 
#         V_t = self.n * self.k * T_c / self.q # Thermal voltage
#         I_0 = self.I_0_ref * (T_c/self.T_c_ref)**3 * \
#             np.exp((self.q * self.E_g) / (self.n * self.k) * (1/self.T_c_ref - 1/T_c))
#         I_ph = (self.I_sc_ref + self.K_i * (T_c - self.T_c_ref)) * (G / self.G_ref)
#         I_d = I_0 * (np.exp(((V/self.N_s) + (I * self.R_s)/self.N_p) / V_t) - 1)
#         I_sh = ((V/self.N_s) + (I * self.R_s)/self.N_p) / self.R_sh 

#         return I - self.N_p * (I_ph - I_d - I_sh)
    
#     def df1_dI(self, I, V, T_c):
#         V_t = self.n * self.k * T_c / self.q # Thermal voltage
#         I_0 = self.I_0_ref * (T_c/self.T_c_ref)**3 * \
#             np.exp((self.q * self.E_g) / (self.n * self.k) * (1/self.T_c_ref - 1/T_c))
        
#         return 1 + (self.N_p * I_0 * np.exp(((V/self.N_s) + (I * self.R_s)/self.N_p) / V_t) * (self.R_s / (V_t * self.N_p))) + (self.R_s / self.R_sh)


#     # Newton-Raphson method to solve for current I for given V, G and T_c
#     def solve_I(self, V, G, T_c):
#         I_ph = (self.I_sc_ref + self.K_i * (T_c - self.T_c_ref)) * (G / self.G_ref)
#         I_guess = self.N_p * I_ph  # Initial guess (Maximum possible current is Iph)
#         tolerance = 1e-6
#         max_iter = 100

#         for i in range(max_iter):
#             f_val = self.f1(I_guess, V, G, T_c)
#             df_dI_val = self.df1_dI(I_guess, V, T_c)

#             if df_dI_val == 0:  # Avoid division by zero
#                 break

#             I_new = I_guess - f_val / df_dI_val

#             if abs(I_new - I_guess) < tolerance:
#                 return max(0, I_new)  # Current cannot be negative

#             I_guess = I_new
        
#         return max(0, I_guess)  # Return the last guess if max iterations reached   
    
#     def plot_IV_PV(self, G, T_c):
#         voltages = np.linspace(0, self.N_s * self.V_oc_ref, 100)
#         currents = [self.solve_I(V, G, T_c) for V in voltages]
#         currents_low_G = [self.solve_I(V, G - 250, T_c) for V in voltages]
#         currents_high_G = [self.solve_I(V, G + 250, T_c) for V in voltages]
#         currents_low_T = [self.solve_I(V, G, T_c - 10) for V in voltages]
#         currents_high_T = [self.solve_I(V, G, T_c + 10) for V in voltages]
#         powers = voltages * currents
#         powers_low_G = voltages * currents_low_G
#         powers_high_G = voltages * currents_high_G
#         powers_low_T = voltages * currents_low_T
#         powers_high_T = voltages * currents_high_T
#         fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
#         ax1.plot(voltages, currents, label='At Standard Conditions')
#         ax1.plot(voltages, currents_low_G, label='Lower Irradiance (G-250 W/m^2)')
#         ax1.plot(voltages, currents_high_G, label='Higher Irradiance (G+250 W/m^2)')
#         ax1.plot(voltages, currents_low_T, label='Lower Temperature (T-10 °C)')
#         ax1.plot(voltages, currents_high_T, label='Higher Temperature (T+10 °C)')
#         ax1.legend()
#         ax1.set_xlabel('Voltage (V)')
#         ax1.set_ylabel('Current (A)')
#         ax1.set_title(f'I-V Curve at G={G} W/m^2 and T_c={T_c-273.15}°C')
#         ax1.grid()
#         ax2.plot(voltages, powers, label='At Standard Conditions')
#         ax2.plot(voltages, powers_low_G, label='Lower Irradiance (G-250 W/m^2)')
#         ax2.plot(voltages, powers_high_G, label='Higher Irradiance (G+250 W/m^2)')
#         ax2.plot(voltages, powers_low_T, label='Lower Temperature (T-10 °C)')
#         ax2.plot(voltages, powers_high_T, label='Higher Temperature (T+10 °C)')
#         ax2.legend()
#         ax2.set_xlabel('Voltage (V)')
#         ax2.set_ylabel('Power (W)')
#         ax2.set_title(f'P-V Curve at G={G} W/m^2 and T_c={T_c-273.15}°C')
#         ax2.grid()
#         plt.show()

# Solar_array = PV_array(N_s=36, N_p=100)
# # Plot I-V and P-V curves at standard conditions    
# Solar_array.plot_IV_PV(Solar_array.G_ref, Solar_array.T_c_ref)




















# PV_array model in shading conditions
class PV_Module:
    def __init__(self, N_cells=36):
        # Physical Constants
        self.q = 1.602176634e-19
        self.k = 1.380649e-23
        self.N_cells = N_cells
        
        # STC Parameters
        self.T_c_ref = 25 + 273.15
        self.G_ref = 1000
        self.E_g = 1.1 * self.q
        self.I_sc_ref = 5.0
        self.K_i = 0.005
        self.V_oc_ref = 0.6  # Per cell
        self.n = 1.3
        self.R_s = 0.01
        
    def calculate_cell_voltage(self, I, G, T_c):
        """Calculates a single cell's voltage for a given current using the single-diode model."""
        # Temperature adjustments
        I_ph = (self.I_sc_ref + self.K_i * (T_c - self.T_c_ref)) * (G / self.G_ref)
        
        # Saturation current
        V_oc_temp = self.V_oc_ref + (-0.002 * (T_c - self.T_c_ref)) # Typical -2mV/K
        I_rs = self.I_sc_ref / (np.exp((self.q * V_oc_temp) / (self.n * self.k * T_c)) - 1)
        I_o = I_rs * (T_c / self.T_c_ref)**3 * np.exp((self.E_g / self.k) * (1/self.T_c_ref - 1/T_c))
        
        # Thermal voltage
        V_t = (self.k * T_c) / self.q
        
        # Avoid log of negative or zero values during implicit solving
        arg = (I_ph - I + I_o) / I_o
        if arg <= 0:
            return 0.0
            
        V_cell = self.n * V_t * np.log(arg) - I * self.R_s
        return max(0.0, V_cell)

    def get_module_voltage(self, I, G, T_c, has_bypass_diode=True):
        """Calculates module voltage including bypass diode behavior."""
        if I > (self.I_sc_ref * (G / self.G_ref)):
            # If current demands exceed short-circuit current, module goes into reverse bias
            if has_bypass_diode:
                return -0.7  # Bypass diode turns on (voltage drop)
            return 0.0
        
        V_cells_total = self.calculate_cell_voltage(I, G, T_c) * self.N_cells
        return V_cells_total

class PV_String:
    def __init__(self, modules):
        self.modules = modules # List of PV_Module instances
        
    def get_string_curves(self, shading_profile, T_c=298.15):
        """
        shading_profile: List of irradiance values [G1, G2, ...] for each module
        """
        # Determine maximum possible current (from the unshaded module)
        max_I = max(mod.I_sc_ref * (G / mod.G_ref) for mod, G in zip(self.modules, shading_profile))
        
        # Generate current array from 0 to max short-circuit current
        current_array = np.linspace(0, max_I * 1.05, 500)
        voltage_array = []
        power_array = []
        
        for I in current_array:
            V_string = 0
            for mod, G in zip(self.modules, shading_profile):
                V_string += mod.get_module_voltage(I, G, T_c)
            
            # Clamp string voltage to 0 (cannot generate negative power overall)
            V_string = max(0.0, V_string)
            voltage_array.append(V_string)
            power_array.append(V_string * I)
            
        return np.array(voltage_array), current_array, np.array(power_array)

# --- Simulation and Plotting ---

# Create a string of 3 modules
modules_string = PV_String([PV_Module(), PV_Module(), PV_Module()])

# Define different sunlight scenarios (Irradiance per module in W/m^2)
scenarios = {
    "Uniform Sunlight (STC)": [1000, 1000, 1000],
    "Mild Partial Shading": [1000, 700, 1000],
    "Severe Partial Shading": [1000, 500, 300]
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

for label, profile in scenarios.items():
    V, I, P = modules_string.get_string_curves(profile, T_c=25+273.15)
    
    # Filter out points where voltage drops to zero to clean up plots
    valid = V > 0
    
    ax1.plot(V[valid], I[valid], label=label, linewidth=2)
    ax2.plot(V[valid], P[valid], label=label, linewidth=2)

ax1.set_xlabel('String Voltage (V)')
ax1.set_ylabel('String Current (A)')
ax1.set_title('I-V Characteristics under Variable Shading')
ax1.grid(True)
ax1.legend()

ax2.set_xlabel('String Voltage (V)')
ax2.set_ylabel('String Power (W)')
ax2.set_title('P-V Characteristics (Multi-Peak Distortion)')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()
    











