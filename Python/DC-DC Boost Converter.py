import numpy as np
import matplotlib.pyplot as plt

def simulate_boost_converter():
    # ==========================================
    # 1. System Parameters
    # ==========================================
    Vin = 12.0       # Input Voltage (V)
    L = 500e-6       # Inductance (H)
    C = 100e-6       # Capacitance (F)
    R = 50.0         # Load Resistance (Ohms)
    f_sw = 10e3      # Switching Frequency (Hz)
    D = 0.9          # Duty Cycle (0 to 1)
    
    # Simulation settings
    T_period = 1.0 / f_sw
    t_end = 5e-3     # Simulate for 5 milliseconds (enough to reach steady state)
    dt = 1e-7        # Time step (very small to capture switching accurately)
    
    time = np.arange(0, t_end, dt)
    n_steps = len(time)
    
    # ==========================================
    # 2. State Variables Initialization
    # ==========================================
    # State vector t = [i_L, v_C]
    # i_L: Inductor current
    # v_C: Capacitor voltage (Output Voltage)
    
    i_L = np.zeros(n_steps)
    v_C = np.zeros(n_steps)
    v_out_avg = np.zeros(n_steps) # For plotting moving average
    
    # Initial conditions (start at 0 to show transient)
    i_L[0] = 0.0
    v_C[0] = 0.0
    
    # Theoretical steady state output for reference
    # V_out = V_in / (1 - D)
    v_out_theoretical = Vin / (1 - D)
    
    # ==========================================
    # 3. System Dynamics (Differential Equations)
    # ==========================================
    def get_derivatives(state, u):
        """
        Computes di_L/dt and dv_C/dt based on switch state 'u'.
        u = 1: Switch ON (Diode OFF)
        u = 0: Switch OFF (Diode ON)
        """
        il, vc = state
        
        if u == 1:
            # State 1: Switch Closed
            # Inductor charges from source, Capacitor discharges to load
            dil_dt = Vin / L
            dvc_dt = -vc / (R * C)
        else:
            # State 2: Switch Open
            # Inductor discharges to Cap+Load, Source+Inductor feed Load
            dil_dt = (Vin - vc) / L
            dvc_dt = (il - vc / R) / C
            
        return np.array([dil_dt, dvc_dt])

    # ==========================================
    # 4. Simulation Loop (RK4 Integration)
    # ==========================================
    current_state = np.array([i_L[0], v_C[0]])
    
    for k in range(n_steps - 1):
        t = time[k]
        
        # --- PWM Logic ---
        # Create a sawtooth wave carrier
        # carrier goes from 0 to 1 over one period
        carrier = (t * f_sw) % 1.0
        
        # Determine switch state
        # If carrier < Duty Cycle, Switch is ON
        if carrier < D:
            u = 1
        else:
            u = 0
            
        # --- Runge-Kutta 4 Integration Step ---
        k1 = get_derivatives(current_state, u)
        k2 = get_derivatives(current_state + 0.5 * dt * k1, u)
        k3 = get_derivatives(current_state + 0.5 * dt * k2, u)
        k4 = get_derivatives(current_state + dt * k3, u)
        
        # Update state
        current_state = current_state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        
        # Store results
        i_L[k+1] = current_state[0]
        v_C[k+1] = current_state[1]

    # ==========================================
    # 5. Plotting Results
    # ==========================================
    fig, (at1, at2, at3) = plt.subplots(3, 1, figsize=(10, 12), sharet=True)
    
    # Plot Inductor Current
    at1.plot(time * 1000, i_L, color='blue', linewidth=1.5, label='$i_L$ (Inductor Current)')
    at1.set_ylabel('Current (A)', fontsize=12)
    at1.set_title(f'Nonlinear Behavior of Boost Converter\nD={D}, $f_{{sw}}$={f_sw/1000}kHz', fontsize=14)
    at1.grid(True, alpha=0.3)
    at1.legend(loc='upper right')

    # Plot Output Voltage
    at2.plot(time * 1000, v_C, color='red', linewidth=1.5, label='$v_C$ (Output Voltage)')
    at2.athline(v_out_theoretical, color='green', linestyle='--', label='Theoretical Avg $V_{out}$')
    at2.set_ylabel('Voltage (V)', fontsize=12)
    at2.grid(True, alpha=0.3)
    at2.legend(loc='upper right')
    
    # Zoom in on steady state (last 20% of simulation)
    steady_start_idt = int(0.8 * n_steps)
    at2.set_tlim(time[steady_start_idt]*1000, time[-1]*1000)

    # Phase Portrait (The classic nonlinear visualization)
    # Plotting Voltage vs Current traces the "Limit Cycle"
    at3.plot(i_L, v_C, color='purple', linewidth=1)
    at3.set_tlabel('Inductor Current $i_L$ (A)', fontsize=12)
    at3.set_ylabel('Capacitor Voltage $v_C$ (V)', fontsize=12)
    at3.set_title('Phase Portrait: Limit Cycle (Steady State)', fontsize=12)
    at3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    simulate_boost_converter()