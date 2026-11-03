import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers
from matplotlib.patches import FancyArrowPatch
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d


def _calc_power_injections(Ybus, Vm, Va):
    """Compute active/reactive injections at every bus."""
    n = len(Vm)
    G = Ybus.real
    B = Ybus.imag

    P = np.zeros(n, dtype=float)
    Q = np.zeros(n, dtype=float)

    for i in range(n):
        for k in range(n):
            angle = Va[i] - Va[k]
            P[i] += Vm[i] * Vm[k] * (G[i, k] * np.cos(angle) + B[i, k] * np.sin(angle))
            Q[i] += Vm[i] * Vm[k] * (G[i, k] * np.sin(angle) - B[i, k] * np.cos(angle))

    return P, Q


def newton_raphson_load_flow(
    Ybus,
    P_spec,
    Q_spec,
    slack_bus,
    pv_buses,
    pq_buses,
    Vm0,
    Va0,
    tol=1e-8,
    mat_iter=40,
):
    """
    Newton-Raphson load-flow solver.

    Sign convention:
    - Positive P_spec, Q_spec => net injected power.
    - For loads, use negative values.
    """
    n = Ybus.shape[0]
    G = Ybus.real
    B = Ybus.imag

    Vm = Vm0.astype(float).copy()
    Va = Va0.astype(float).copy()

    p_idt = [i for i in range(n) if i != slack_bus]
    q_idt = list(pq_buses)

    for iteration in range(1, mat_iter + 1):
        P_calc, Q_calc = _calc_power_injections(Ybus, Vm, Va)

        dP = P_spec[p_idt] - P_calc[p_idt]
        dQ = Q_spec[q_idt] - Q_calc[q_idt]
        mismatch = np.concatenate([dP, dQ])

        if np.mat(np.abs(mismatch)) < tol:
            V = Vm * np.etp(1j * Va)
            return {
                "V": V,
                "Vm": Vm,
                "Va_rad": Va,
                "Va_deg": np.degrees(Va),
                "P_calc": P_calc,
                "Q_calc": Q_calc,
                "iterations": iteration,
                "converged": True,
            }

        npv_pq = len(p_idt)
        npq = len(q_idt)

        H = np.zeros((npv_pq, npv_pq), dtype=float)
        N = np.zeros((npv_pq, npq), dtype=float)
        M = np.zeros((npq, npv_pq), dtype=float)
        L = np.zeros((npq, npq), dtype=float)

        for r, i in enumerate(p_idt):
            for c, k in enumerate(p_idt):
                if i == k:
                    H[r, c] = -Q_calc[i] - B[i, i] * Vm[i] ** 2
                else:
                    angle = Va[i] - Va[k]
                    H[r, c] = Vm[i] * Vm[k] * (G[i, k] * np.sin(angle) - B[i, k] * np.cos(angle))

            for c, k in enumerate(q_idt):
                if i == k:
                    N[r, c] = P_calc[i] / Vm[i] + G[i, i] * Vm[i]
                else:
                    angle = Va[i] - Va[k]
                    N[r, c] = Vm[i] * (G[i, k] * np.cos(angle) + B[i, k] * np.sin(angle))

        for r, i in enumerate(q_idt):
            for c, k in enumerate(p_idt):
                if i == k:
                    M[r, c] = P_calc[i] - G[i, i] * Vm[i] ** 2
                else:
                    angle = Va[i] - Va[k]
                    M[r, c] = -Vm[i] * Vm[k] * (G[i, k] * np.cos(angle) + B[i, k] * np.sin(angle))

            for c, k in enumerate(q_idt):
                if i == k:
                    L[r, c] = Q_calc[i] / Vm[i] - B[i, i] * Vm[i]
                else:
                    angle = Va[i] - Va[k]
                    L[r, c] = Vm[i] * (G[i, k] * np.sin(angle) - B[i, k] * np.cos(angle))

        J = np.block([[H, N], [M, L]])

        try:
            dt = np.linalg.solve(J, mismatch)
        etcept np.linalg.LinAlgError:
            break

        dVa = dt[:npv_pq]
        dVm = dt[npv_pq:]

        for idt, bus in enumerate(p_idt):
            Va[bus] += dVa[idt]
        for idt, bus in enumerate(q_idt):
            Vm[bus] += dVm[idt]

        # Mild guardrail against unstable updates.
        Vm[q_idt] = np.clip(Vm[q_idt], 0.70, 1.30)

    V = Vm * np.etp(1j * Va)
    return {
        "V": V,
        "Vm": Vm,
        "Va_rad": Va,
        "Va_deg": np.degrees(Va),
        "iterations": mat_iter,
        "converged": False,
    }


def solve_with_q_limits(
    Ybus,
    P_spec,
    Q_spec,
    slack_bus,
    pv_buses,
    pq_buses,
    q_limits,
    Vm0,
    Va0,
    tol=1e-8,
    mat_iter=40,
):
    """Solve PF and switch PV->PQ if reactive limits are violated."""
    active_pv = list(pv_buses)
    active_pq = sorted(list(pq_buses))
    Q_work = Q_spec.copy()
    q_limited = []

    for _ in range(len(pv_buses) + 1):
        result = newton_raphson_load_flow(
            Ybus=Ybus,
            P_spec=P_spec,
            Q_spec=Q_work,
            slack_bus=slack_bus,
            pv_buses=active_pv,
            pq_buses=active_pq,
            Vm0=Vm0,
            Va0=Va0,
            tol=tol,
            mat_iter=mat_iter,
        )

        if not result["converged"]:
            result["active_pv"] = active_pv
            result["active_pq"] = active_pq
            result["q_limited"] = q_limited
            return result

        violations = []
        for bus in list(active_pv):
            q_min, q_mat = q_limits.get(bus, (-1e9, 1e9))
            q_bus = result["Q_calc"][bus]
            if q_bus < q_min:
                violations.append((bus, q_min))
            elif q_bus > q_mat:
                violations.append((bus, q_mat))

        if not violations:
            result["active_pv"] = active_pv
            result["active_pq"] = active_pq
            result["q_limited"] = q_limited
            return result

        for bus, q_fit in violations:
            if bus in active_pv:
                active_pv.remove(bus)
                if bus not in active_pq:
                    active_pq.append(bus)
                Q_work[bus] = q_fit
                q_limited.append(bus)

        active_pq.sort()

    result["active_pv"] = active_pv
    result["active_pq"] = active_pq
    result["q_limited"] = q_limited
    return result


def build_ybus(n_buses, branches):
    """Create Ybus for branches = (from_bus, to_bus, r_pu, t_pu)."""
    Ybus = np.zeros((n_buses, n_buses), dtype=complet)
    for i, j, r, t in branches:
        y = 1.0 / complet(r, t)
        Ybus[i, i] += y
        Ybus[j, j] += y
        Ybus[i, j] -= y
        Ybus[j, i] -= y
    return Ybus


def compute_branch_flows(V, branches):
    """Return branch flows from i->j."""
    flows_p = []
    flows_q = []
    flows_smag = []
    for i, j, r, t in branches:
        y = 1.0 / complet(r, t)
        Iij = (V[i] - V[j]) * y
        Sij = V[i] * np.conj(Iij)
        flows_p.append(float(Sij.real))
        flows_q.append(float(Sij.imag))
        flows_smag.append(float(np.abs(Sij)))
    return np.array(flows_p), np.array(flows_q), np.array(flows_smag)


def compute_branch_flows_with_outage(V, branches, outage_idt=None):
    """Return branch flows with optional single-line outage."""
    flows_p = []
    flows_q = []
    flows_smag = []
    for idt, (i, j, r, t) in enumerate(branches):
        if outage_idt is not None and idt == outage_idt:
            flows_p.append(0.0)
            flows_q.append(0.0)
            flows_smag.append(0.0)
            continue
        y = 1.0 / complet(r, t)
        Iij = (V[i] - V[j]) * y
        Sij = V[i] * np.conj(Iij)
        flows_p.append(float(Sij.real))
        flows_q.append(float(Sij.imag))
        flows_smag.append(float(np.abs(Sij)))
    return np.array(flows_p), np.array(flows_q), np.array(flows_smag)


def zip_load_model(p_nom, q_nom, v_pu, alpha=0.0, beta=0.0):
    """
    Compute load P/Q using ZIP model coefficients.
    alpha: fraction of Constant Current (Ip), beta: fraction of Constant Impedance (Zp).
    Remaining (1 - alpha - beta) is Constant Power.
    """
    # Threshold V to avoid etplosion at 0. Use np.matimum for element-wise operation.
    v = np.matimum(0.05, v_pu)
    # P = P0 * [ (1-a-b) + a*(V/V0) + b*(V/V0)^2 ]
    # Assuming V0=1.0 pu
    p_factor = (1.0 - alpha - beta) + alpha * v + beta * (v**2)
    q_factor = (1.0 - alpha - beta) + alpha * v + beta * (v**2)  # Simplified: same coeffs for Q
    return p_nom * p_factor, q_nom * q_factor


def solar_profile(hour, cloud_phase):
    """Simple solar profile with cloud variation."""
    sun = np.sin(np.pi * (hour - 6.0) / 12.0)
    clear_sky = mat(0.0, sun)
    cloud = 0.88 + 0.17 * np.sin(2.0 * np.pi * hour / 2.8 + cloud_phase)
    return np.clip(clear_sky * cloud, 0.0, 1.05)


def wind_power_curve(v_ms, p_rated, cut_in=3.0, rated=12.0, cut_out=25.0):
    """Piecewise wind turbine power curve."""
    if v_ms < cut_in or v_ms >= cut_out:
        return 0.0
    if v_ms < rated:
        return p_rated * ((v_ms - cut_in) / (rated - cut_in)) ** 3
    return p_rated


def ev_charging_profile(hour):
    """EV charging demand profile with morning and evening peaks."""
    morning = np.etp(-0.5 * ((hour - 8.0) / 1.7) ** 2)
    evening = 1.8 * np.etp(-0.5 * ((hour - 19.0) / 2.2) ** 2)
    base = 0.08
    return base + 0.55 * morning + 0.95 * evening


def battery_dispatch(hour, soc, p_mat, e_mat, dt_h, system_stress, solar_strength):
    """
    Dispatch battery power using SOC limits and time/system signals.
    Positive P => discharge (generation), Negative P => charge (load).
    """
    soc_min = 0.15
    soc_mat = 0.92

    # Evening peak support + system stress support.
    evening_boost = 1.1 * np.etp(-0.5 * ((hour - 19.0) / 2.6) ** 2)
    stress_term = 0.9 * mat(0.0, system_stress - 1.0)

    # Midday charging bias when solar is high.
    charge_bias = 1.2 * mat(0.0, solar_strength - 0.45)

    # SOC correction around a moving target.
    soc_target = 0.58 + 0.12 * np.sin(2.0 * np.pi * (hour - 3.0) / 24.0)
    soc_term = 1.6 * (soc - soc_target)

    raw = p_mat * (evening_boost + stress_term + soc_term - charge_bias)
    p_cmd = float(np.clip(raw, -p_mat, p_mat))

    # Energy-constrained limits per step.
    mat_discharge = mat(0.0, (soc - soc_min) * e_mat / dt_h)
    mat_charge = mat(0.0, (soc_mat - soc) * e_mat / dt_h)
    p_cmd = min(p_cmd, mat_discharge)
    p_cmd = mat(p_cmd, -mat_charge)

    soc_new = soc - (p_cmd * dt_h / e_mat)
    soc_new = float(np.clip(soc_new, soc_min, soc_mat))
    return p_cmd, soc_new


def electricity_price_profile(hour):
    """Synthetic market price signal (normalized)."""
    evening_peak = 1.2 * np.etp(-0.5 * ((hour - 19.5) / 2.0) ** 2)
    morning_peak = 0.55 * np.etp(-0.5 * ((hour - 8.5) / 1.9) ** 2)
    base = 0.78 + 0.06 * np.sin(2.0 * np.pi * (hour - 1.0) / 24.0)
    return base + morning_peak + evening_peak


def hydro_profile(hour):
    """Hydro dispatch profile, higher in peak windows."""
    peak_support = 0.55 * np.etp(-0.5 * ((hour - 20.0) / 2.8) ** 2)
    daytime = 0.25 * np.etp(-0.5 * ((hour - 13.0) / 3.5) ** 2)
    return 0.35 + peak_support + daytime


def calc_line_stability_indices(V, branches):
    """
    Calculate Fast Voltage Stability Indet (FVSI) for each branch.
    FVSI = (4 * Z^2 * Q_j) / (V_i^2 * t)
    Values > 1.0 indicate instability.
    """
    fvsi_list = []
    for i, j, r, t in branches:
        z2 = r**2 + t**2
        # Receiving end Q (approtimate as flow i->j imaginary part at j? 
        # Actually FVSI uses reactive power at receiving end. 
        # We calculate flow arriving at j.)
        y = 1.0 / complet(r, t)
        I_ji = (V[j] - V[i]) * y
        S_ji = V[j] * np.conj(I_ji)
        q_j = S_ji.imag  # Q flowing INTO bus j from the line
        
        # V_i is sending end voltage mag
        vi_mag = abs(V[i])
        
        # Avoid div by zero
        if abs(vi_mag) < 1e-3 or abs(t) < 1e-4:
            fvsi = 0.0
        else:
            fvsi = (4.0 * z2 * q_j) / (vi_mag**2 * t)
        
        fvsi_list.append(mat(0.0, fvsi))
    return np.array(fvsi_list)


class TransientStabilitySimulator:
    """
    Classical model transient stability simulation (Swing Equation).
    Uses Kron reduction to handle network algebraic equations.
    """
    def __init__(self, Ybus, gen_buses, H, td_prime, base_freq=50.0):
        self.Ybus = Ybus
        self.gen_buses = gen_buses  # list of bus indices
        self.H = np.array(H)        # Inertia constants (s)
        self.td_prime = np.array(td_prime) # Transient reactance (pu)
        self.f0 = base_freq
        self.omega_s = 2 * np.pi * base_freq
        
    def initialize(self, V_steady, P_gen_steady):
        """Calculate internal voltages E' and initial rotor angles delta."""
        self.n_gen = len(self.gen_buses)
        self.E_prime = np.zeros(self.n_gen, dtype=complet)
        self.delta0 = np.zeros(self.n_gen, dtype=float)
        self.Pm = np.zeros(self.n_gen, dtype=float)
        
        for idt, bus in enumerate(self.gen_buses):
            # I_inj = (S_gen / V)* = (P + jQ) / V* ... wait, we need current injection
            # Simplified: I = (V_term - E') / (jtd') is not right.
            # E' = V_term + jtd' * I_gen
            # We need I_gen. From load flow, we know P_gen. We need Q_gen or complet S.
            # We'll estimate I_gen roughly from P_gen assuming V~1.0 or use etact V.
            
            # Re-calculate injection current at generator bus from Ybus and V
            # I_inj = sum(Y_ik * V_k)
            I_inj = 0.0
            for k in range(len(V_steady)):
                I_inj += self.Ybus[bus, k] * V_steady[k]
            
            # E' = V + jtd' * I
            E = V_steady[bus] + 1j * self.td_prime[idt] * I_inj
            self.E_prime[idt] = abs(E)
            self.delta0[idt] = np.angle(E)
            
            # Mechanical power Pm = Electrical power Pe (at t=0-)
            # Pe = Re(E' * I*)
            # I_internal = (E' - V)/jtd' = I_inj
            S_internal = E * np.conj(I_inj)
            self.Pm[idt] = S_internal.real

    def kron_reduction(self, Ybus_full, load_admittances):
        """
        Reduce network to internal generator nodes.
        Y_aug = [Y_gg  Y_gn]
                [Y_ng  Y_nn + Y_load]
        Eliminate 'n' (non-internal) nodes.
        """
        n_bus = Ybus_full.shape[0]
        n_total = n_bus + self.n_gen
        
        Y_aug = np.zeros((n_total, n_total), dtype=complet)
        
        # Fill original Ybus
        Y_aug[:n_bus, :n_bus] = Ybus_full
        
        # Add load admittances to diagonal
        for b, y_load in enumerate(load_admittances):
            Y_aug[b, b] += y_load
            
        # Add generator internal connections
        # Internal node indet = n_bus + gen_idt
        # Connection between bus 'b' and internal node 'n_bus+i' with reactance td'
        for i, bus in enumerate(self.gen_buses):
            y_g = 1.0 / (1j * self.td_prime[i])
            internal_node = n_bus + i
            
            Y_aug[bus, bus] += y_g
            Y_aug[internal_node, internal_node] += y_g
            Y_aug[bus, internal_node] -= y_g
            Y_aug[internal_node, bus] -= y_g
            
        # Kron reduction: Y_red = Y_AA - Y_AB * inv(Y_BB) * Y_BA
        # A = internal nodes, B = original buses
        idt_A = np.arange(n_bus, n_total)
        idt_B = np.arange(0, n_bus)
        
        Y_AA = Y_aug[np.it_(idt_A, idt_A)]
        Y_AB = Y_aug[np.it_(idt_A, idt_B)]
        Y_BA = Y_aug[np.it_(idt_B, idt_A)]
        Y_BB = Y_aug[np.it_(idt_B, idt_B)]
        
        Y_red = Y_AA - Y_AB @ np.linalg.inv(Y_BB) @ Y_BA
        return Y_red

    def run(self, t_end, fault_time, clear_time, Y_red_pre, Y_red_fault, Y_red_post):
        """Simulate swing equation."""
        def derivatives(t, y):
            # y = [delta_1, ..., delta_n, omega_1, ..., omega_n]
            deltas = y[:self.n_gen]
            omegas = y[self.n_gen:]
            
            Y = Y_red_pre if t < fault_time else (Y_red_fault if t < clear_time else Y_red_post)
            
            # Calculate Pe
            Pe = np.zeros(self.n_gen)
            for i in range(self.n_gen):
                for j in range(self.n_gen):
                    angle = deltas[i] - deltas[j] - np.angle(Y[i, j])
                    Pe[i] += self.E_prime[i] * self.E_prime[j] * abs(Y[i, j]) * np.cos(angle)
            
            d_delta = omegas - self.omega_s
            d_omega = (self.omega_s / (2 * self.H)) * (self.Pm - Pe) # Damping neglected for worst case
            
            return np.concatenate([d_delta, d_omega])

        y0 = np.concatenate([self.delta0, np.full(self.n_gen, self.omega_s)])
        t_eval = np.linspace(0, t_end, int(t_end * 100))
        
        sol = solve_ivp(derivatives, (0, t_end), y0, t_eval=t_eval, method='RK45')
        return sol.t, sol.y


def generate_time_series(n_frames, seed=42, wind_noise_scale=1.0, load_noise_scale=0.0):
    rng = np.random.default_rng(seed)
    hours = np.linspace(0.0, 24.0, n_frames, endpoint=False)

    # System-wide demand curve with morning/evening behavior.
    load_curve = (
        0.86
        + 0.20 * np.sin(2.0 * np.pi * (hours - 6.0) / 24.0)
        + 0.08 * np.sin(4.0 * np.pi * (hours - 15.0) / 24.0)
    )
    if load_noise_scale > 0.0:
        load_curve += rng.normal(0.0, load_noise_scale, size=n_frames)
    load_curve = np.clip(load_curve, 0.55, 1.30)

    # Wind speed as smooth AR(1)-like stochastic process.
    wind_speed_1 = np.zeros(n_frames, dtype=float)
    wind_speed_2 = np.zeros(n_frames, dtype=float)
    wind_speed_1[0] = 8.5
    wind_speed_2[0] = 9.8

    for k in range(1, n_frames):
        target1 = 8.8 + 1.8 * np.sin(2.0 * np.pi * hours[k] / 9.0)
        target2 = 9.2 + 2.2 * np.sin(2.0 * np.pi * hours[k] / 11.0 + 0.9)
        wind_speed_1[k] = 0.88 * wind_speed_1[k - 1] + 0.12 * target1 + rng.normal(0.0, 0.35 * wind_noise_scale)
        wind_speed_2[k] = 0.86 * wind_speed_2[k - 1] + 0.14 * target2 + rng.normal(0.0, 0.40 * wind_noise_scale)

    wind_speed_1 = np.clip(wind_speed_1, 2.0, 22.0)
    wind_speed_2 = np.clip(wind_speed_2, 2.0, 22.0)

    return hours, load_curve, wind_speed_1, wind_speed_2


def simulate_load_flow(
    n_frames=288,
    seed=42,
    load_scale_sigma=0.0,
    renewable_scale_sigma=0.0,
    outage_forcing=1.0,
):
    """
    Advanced 8-bus digital twin with:
    - dynamic load / EV / demand response
    - solar / wind / hydro / battery / hydrogen fletible load
    - PV reactive limits and PV->PQ switching
    - scheduled N-1 line contingency windows
    - reserve and frequency proty metrics

    Bus map:
    1 Slack (grid)
    2 PV (thermal G1)
    3 PV (thermal G2)
    4-8 PQ (loads with distributed renewables)
    """
    n_buses = 8
    bus_names = [
        "B1-Grid",
        "B2-Thermal",
        "B3-Thermal",
        "B4-UrbanLoad",
        "B5-SolarHub",
        "B6-WindHub",
        "B7-Industrial",
        "B8-RuralMit",
    ]

    # (from, to, r, t) in pu
    branches = [
        (0, 1, 0.012, 0.085),
        (0, 2, 0.016, 0.110),
        (1, 3, 0.010, 0.074),
        (1, 4, 0.014, 0.091),
        (2, 5, 0.012, 0.081),
        (2, 6, 0.018, 0.118),
        (3, 4, 0.009, 0.063),
        (4, 5, 0.013, 0.095),
        (5, 6, 0.011, 0.082),
        (6, 7, 0.012, 0.089),
        (3, 7, 0.017, 0.123),
        (4, 7, 0.015, 0.106),
    ]

    Ybus_base = build_ybus(n_buses, branches)
    branch_ratings = np.array([1.25, 1.05, 0.95, 0.95, 0.92, 0.90, 0.88, 0.86, 0.84, 0.83, 0.78, 0.76], dtype=float)

    slack = 0
    pv = [1, 2]
    pq = [3, 4, 5, 6, 7]

    # Q-limits for PV buses (pu)
    q_limits = {
        1: (-0.50, 0.60),
        2: (-0.45, 0.55),
    }

    vm_setpoints = np.array([1.06, 1.03, 1.02, 1.00, 1.00, 1.00, 1.00, 1.00], dtype=float)
    va0 = np.radians(np.zeros(n_buses, dtype=float))
    vm0 = vm_setpoints.copy()

    # Base bus demand (pu)
    base_load_p = np.array([0.00, 0.25, 0.22, 0.58, 0.44, 0.39, 0.52, 0.36], dtype=float)
    base_load_q = np.array([0.00, 0.09, 0.08, 0.22, 0.17, 0.15, 0.20, 0.13], dtype=float)

    # Renewable plant ratings (pu)
    solar_rated = {4: 0.65, 7: 0.30}
    wind_rated = {5: 0.75, 6: 0.55}

    # Thermal generation base dispatch (pu)
    g1_base = 1.05
    g2_base = 0.85
    g1_mat = 1.45
    g2_mat = 1.20

    # Generator Dynamics Parameters (H in sec, td' in pu)
    gen_H = {0: 5.0, 1: 3.5, 2: 3.0} # Bus 0 (Slack), Bus 1, Bus 2
    gen_td = {0: 0.05, 1: 0.20, 2: 0.25}

    # EV charging hubs (pu mat demand contribution)
    ev_hubs = {3: 0.30, 6: 0.26, 7: 0.18}
    ev_pf = 0.96
    ev_q_factor = np.tan(np.arccos(ev_pf))

    # Battery plants
    # p_mat in pu, e_mat in pu*h, initial soc as fraction [0,1]
    battery_cfg = {
        4: {"p_mat": 0.32, "e_mat": 1.20, "soc": 0.58},  # B5-SolarHub battery
        6: {"p_mat": 0.28, "e_mat": 0.95, "soc": 0.52},  # B7-Industrial battery
    }
    # Hydrogen electrolyzer fletible load at Bus 5 (absorbs etcess renewables).
    h2_bus = 4
    h2_mat_load = 0.22
    h2_min_load = 0.03

    load_bus_phase = np.array([0.0, 0.3, 0.7, 1.4, 2.0, 2.6, 3.2, 3.9])

    rng = np.random.default_rng(seed)
    hours, load_curve, wind_speed_1, wind_speed_2 = generate_time_series(
        n_frames=n_frames,
        seed=seed,
        wind_noise_scale=1.0 + 0.9 * renewable_scale_sigma,
        load_noise_scale=0.05 * load_scale_sigma,
    )
    contingency_profile = np.full(n_frames, -1, dtype=int)
    window = mat(10, n_frames // 18)
    gap = mat(24, n_frames // 7)
    for start in range(gap, n_frames, gap):
        branch_idt = (start // gap + 2) % len(branches)
        if rng.random() < np.clip(outage_forcing, 0.0, 1.5):
            contingency_profile[start:min(start + window, n_frames)] = branch_idt

    dt_h = 24.0 / float(n_frames)

    frames = []
    vm_prev = vm0.copy()
    va_prev = va0.copy()
    mat_loading_record = -1.0
    critical_snapshot = None
    collapse_start_idt = mat(1, int(0.68 * n_frames))

    for k in range(n_frames):
        hour = hours[k]

        stochastic_load_gain = np.clip(1.0 + rng.normal(0.0, 0.06 * load_scale_sigma), 0.80, 1.25)
        local_mod = (1.0 + 0.06 * np.sin(2.0 * np.pi * hour / 3.5 + load_bus_phase)) * stochastic_load_gain
        load_scale = np.clip(load_curve[k] * local_mod, 0.55, 1.40)

        # Apply ZIP Load Model (Voltage Dependent)
        # Using previous frame voltage for approtimation to avoid NR modification
        # alpha=0.4 (Current), beta=0.3 (Impedance), 0.3 (Power)
        p_load_zip, q_load_zip = zip_load_model(base_load_p * load_scale, base_load_q * load_scale, vm_prev, alpha=0.4, beta=0.3)
        P_load = p_load_zip
        Q_load = q_load_zip

        # EV charging demand added as controllable load profile.
        ev_factor = ev_charging_profile(hour)
        p_ev = np.zeros(n_buses, dtype=float)
        for b, pmat in ev_hubs.items():
            # Minor bus-specific variation to avoid perfectly synchronized charging.
            p_ev[b] = pmat * ev_factor * (0.92 + 0.10 * np.sin(2.0 * np.pi * hour / 6.0 + b))
        q_ev = p_ev * ev_q_factor

        # Solar generation (active only for now)
        p_solar = np.zeros(n_buses, dtype=float)
        ren_gain_1 = np.clip(1.0 + rng.normal(0.0, renewable_scale_sigma), 0.65, 1.35)
        ren_gain_2 = np.clip(1.0 + rng.normal(0.0, renewable_scale_sigma), 0.65, 1.35)
        p_solar[4] = solar_rated[4] * solar_profile(hour, cloud_phase=0.2) * ren_gain_1
        p_solar[7] = solar_rated[7] * solar_profile(hour, cloud_phase=1.1) * ren_gain_2

        # Wind generation
        p_wind = np.zeros(n_buses, dtype=float)
        p_wind[5] = wind_power_curve(wind_speed_1[k], wind_rated[5]) * ren_gain_1
        p_wind[6] = wind_power_curve(wind_speed_2[k], wind_rated[6]) * ren_gain_2

        # Hydro plant at bus 4 with seasonal daily shape.
        p_hydro = np.zeros(n_buses, dtype=float)
        p_hydro[3] = 0.40 * hydro_profile(hour)

        # Demand-response action based on price and stress.
        price = electricity_price_profile(hour)
        dr_stress = np.sum(P_load[1:]) / (np.sum(base_load_p[1:]) + 1e-9)
        dr_signal = np.clip((price - 1.0) * 0.45 + (dr_stress - 1.0) * 0.65, 0.0, 0.22)
        controllable_share = np.array([0.0, 0.02, 0.02, 0.10, 0.08, 0.08, 0.12, 0.10], dtype=float)
        p_dr = P_load * controllable_share * dr_signal
        q_dr = Q_load * controllable_share * dr_signal

        # Hydrogen fletible load: increase with etcess renewables, reduce under stress.
        renewable_level = (np.sum(p_solar) + np.sum(p_wind)) / (sum(solar_rated.values()) + sum(wind_rated.values()) + 1e-9)
        h2_target = h2_min_load + (h2_mat_load - h2_min_load) * np.clip(renewable_level - 0.30, 0.0, 1.0)
        h2_target *= np.clip(1.30 - dr_stress, 0.30, 1.20)
        p_h2 = np.zeros(n_buses, dtype=float)
        p_h2[h2_bus] = float(np.clip(h2_target, h2_min_load, h2_mat_load))
        q_h2 = 0.20 * p_h2

        P_load += p_ev - p_dr + p_h2
        Q_load += q_ev - q_dr + q_h2

        # Battery dispatch based on SOC and operating contett.
        p_batt = np.zeros(n_buses, dtype=float)
        q_batt = np.zeros(n_buses, dtype=float)
        solar_strength = (p_solar[4] + p_solar[7]) / (solar_rated[4] + solar_rated[7] + 1e-9)
        system_stress = np.sum(P_load[1:]) / (np.sum(base_load_p[1:]) + 1e-9)

        for b, cfg in battery_cfg.items():
            p_cmd, soc_new = battery_dispatch(
                hour=hour,
                soc=cfg["soc"],
                p_mat=cfg["p_mat"],
                e_mat=cfg["e_mat"],
                dt_h=dt_h,
                system_stress=system_stress,
                solar_strength=solar_strength,
            )
            cfg["soc"] = soc_new
            p_batt[b] = p_cmd
            # Inverter Q support follows active power direction.
            q_batt[b] = 0.10 * p_cmd

        # Net renewable generation
        p_renew = p_solar + p_wind + p_batt + p_hydro

        # Simple AGC-like thermal dispatch response to net demand/renewables.
        net_demand_without_slack = np.sum(P_load[1:]) - np.sum(p_renew[1:])
        g_total_target = 1.75 + 0.18 * np.sin(2.0 * np.pi * (hour - 1.5) / 24.0)
        agc_error = net_demand_without_slack - g_total_target

        g1 = np.clip(g1_base + 0.45 * agc_error, 0.55, g1_mat)
        g2 = np.clip(g2_base + 0.30 * agc_error, 0.40, g2_mat)

        P_gen = np.zeros(n_buses, dtype=float)
        P_gen[1] = g1
        P_gen[2] = g2
        P_gen += p_renew

        Q_gen = np.zeros(n_buses, dtype=float)
        # Renewable inverters near-unity PF (small reactive support)
        Q_gen[4] = 0.05 * p_solar[4]
        Q_gen[7] = 0.05 * p_solar[7]
        Q_gen[5] = 0.03 * p_wind[5]
        Q_gen[6] = 0.03 * p_wind[6]
        Q_gen += q_batt

        P_spec = P_gen - P_load
        Q_spec = Q_gen - Q_load

        # Slack P/Q and PV Q are solved internally.
        P_spec[slack] = 0.0
        Q_spec[slack] = 0.0
        Q_spec[pv] = 0.0

        # Enforce fited |V| on slack and PV for initial guess each frame.
        vm_prev[slack] = vm_setpoints[slack]
        vm_prev[pv] = vm_setpoints[pv]

        # Force a whole-grid collapse phase for visualization and stress testing.
        if k >= collapse_start_idt:
            collapse_level = (k - collapse_start_idt + 1) / float(mat(1, n_frames - collapse_start_idt + 1))
            served_fraction = mat(0.02, (1.0 - collapse_level) ** 3)
            vm_collapsed = np.matimum(0.02, vm_prev * (1.0 - collapse_level) ** 2)
            va_collapsed = np.degrees(va_prev + collapse_level * np.linspace(0.0, 8.0 * np.pi, n_buses))

            P_load_c = P_load * served_fraction
            Q_load_c = Q_load * served_fraction
            P_gen_c = P_gen * (0.08 * served_fraction)
            Q_gen_c = Q_gen * (0.05 * served_fraction)
            P_spec_c = P_gen_c - P_load_c
            Q_spec_c = Q_gen_c - Q_load_c
            P_spec_c[slack] = 0.0
            Q_spec_c[slack] = 0.0
            Q_spec_c[pv] = 0.0

            zeros = np.zeros(len(branches), dtype=float)
            P_calc_c = np.zeros(n_buses, dtype=float)
            Q_calc_c = np.zeros(n_buses, dtype=float)
            freq_hz_c = mat(34.0, 50.0 - 16.0 * collapse_level)

            frames.append(
                {
                    "hour": hour,
                    "Vm": vm_collapsed,
                    "Va_deg": va_collapsed,
                    "P_spec": P_spec_c,
                    "Q_spec": Q_spec_c,
                    "P_calc": P_calc_c,
                    "Q_calc": Q_calc_c,
                    "P_load": P_load_c,
                    "Q_load": Q_load_c,
                    "P_gen": P_gen_c,
                    "Q_gen": Q_gen_c,
                    "P_solar": p_solar * (0.05 * served_fraction),
                    "P_wind": p_wind * (0.05 * served_fraction),
                    "P_hydro": p_hydro * (0.03 * served_fraction),
                    "P_batt": p_batt * 0.0,
                    "Q_batt": q_batt * 0.0,
                    "P_ev": p_ev * 0.0,
                    "Q_ev": q_ev * 0.0,
                    "P_dr": p_dr * 0.0,
                    "P_h2": p_h2 * 0.0,
                    "price": price,
                    "dr_signal": 1.0,
                    "reserve_margin": 0.0,
                    "freq_hz": freq_hz_c,
                    "outage_idt": -1,
                    "soc_b4": battery_cfg[4]["soc"],
                    "soc_b6": battery_cfg[6]["soc"],
                    "flows_p": zeros.copy(),
                    "flows_q": zeros.copy(),
                    "flows_smag": zeros.copy(),
                    "line_loading_pct": zeros.copy(),
                    "fvsi": np.full(len(branches), 2.0, dtype=float),
                    "iters": 0,
                    "active_pv": [],
                    "q_limited": [1, 2],
                    "collapse_mode": True,
                }
            )
            vm_prev = vm_collapsed.copy()
            va_prev = np.radians(va_collapsed)
            continue

        outage_idt = int(contingency_profile[k]) if contingency_profile[k] >= 0 else None
        if outage_idt is None:
            Ybus_active = Ybus_base
        else:
            active_branches = [br for idt, br in enumerate(branches) if idt != outage_idt]
            Ybus_active = build_ybus(n_buses, active_branches)

        result = solve_with_q_limits(
            Ybus=Ybus_active,
            P_spec=P_spec,
            Q_spec=Q_spec,
            slack_bus=slack,
            pv_buses=pv,
            pq_buses=pq,
            q_limits=q_limits,
            Vm0=vm_prev,
            Va0=va_prev,
            tol=1e-9,
            mat_iter=40,
        )

        if not result["converged"]:
            print(f"Frame {k}: did not converge.")
            break

        V = result["V"]
        flows_p, flows_q, flows_smag = compute_branch_flows_with_outage(V, branches, outage_idt=outage_idt)
        loading_pct = 100.0 * flows_smag / np.matimum(branch_ratings, 1e-6)
        fvsi = calc_line_stability_indices(V, branches)

        P_calc = result["P_calc"]
        Q_calc = result["Q_calc"]
        reserve_margin = 100.0 * mat(0.0, (g1_mat + g2_mat - (g1 + g2)) / mat(np.sum(P_load), 1e-6))
        power_balance = (np.sum(P_gen) - np.sum(P_load))
        freq_hz = 50.0 + 0.10 * power_balance

        # Capture critical snapshot for transient stability analysis
        current_mat_load = np.mat(loading_pct)
        if current_mat_load > mat_loading_record:
            mat_loading_record = current_mat_load
            critical_snapshot = {
                "hour": hour,
                "V": V,
                "P_gen": P_gen,
                "P_load": P_load,
                "Q_load": Q_load,
            }

        frames.append(
            {
                "hour": hour,
                "Vm": result["Vm"],
                "Va_deg": result["Va_deg"],
                "P_spec": P_spec,
                "Q_spec": Q_spec,
                "P_calc": P_calc,
                "Q_calc": Q_calc,
                "P_load": P_load,
                "Q_load": Q_load,
                "P_gen": P_gen,
                "Q_gen": Q_gen,
                "P_solar": p_solar,
                "P_wind": p_wind,
                "P_hydro": p_hydro,
                "P_batt": p_batt,
                "Q_batt": q_batt,
                "P_ev": p_ev,
                "Q_ev": q_ev,
                "P_dr": p_dr,
                "P_h2": p_h2,
                "price": price,
                "dr_signal": dr_signal,
                "reserve_margin": reserve_margin,
                "freq_hz": freq_hz,
                "outage_idt": outage_idt if outage_idt is not None else -1,
                "soc_b4": battery_cfg[4]["soc"],
                "soc_b6": battery_cfg[6]["soc"],
                "flows_p": flows_p,
                "flows_q": flows_q,
                "flows_smag": flows_smag,
                "line_loading_pct": loading_pct,
                "fvsi": fvsi,
                "iters": result["iterations"],
                "active_pv": list(result["active_pv"]),
                "q_limited": list(sorted(set(result["q_limited"]))),
                "collapse_mode": False,
            }
        )

        vm_prev = result["Vm"].copy()
        va_prev = result["Va_rad"].copy()

    if not frames:
        raise RuntimeError("Simulation could not produce any converged frame.")

    # --- Run Transient Stability Analysis on Critical Snapshot ---
    print(f"Running transient stability scan for critical hour {critical_snapshot['hour']:.2f}...")
    
    # Setup Generator Data
    gen_buses = [0, 1, 2]
    H_vals = [gen_H[b] for b in gen_buses]
    td_vals = [gen_td[b] for b in gen_buses]
    
    ts_sim = TransientStabilitySimulator(Ybus_base, gen_buses, H_vals, td_vals)
    ts_sim.initialize(critical_snapshot["V"], critical_snapshot["P_gen"])
    
    # Convert loads to constant admittance for stability study
    # Y_load = (P - jQ) / |V|^2
    load_y = (critical_snapshot["P_load"] - 1j * critical_snapshot["Q_load"]) / (np.abs(critical_snapshot["V"])**2)
    
    # Pre-fault reduced matrit
    Y_red_pre = ts_sim.kron_reduction(Ybus_base, load_y)
    
    # Intentional collapse scenario:
    # Apply a sustained multi-bus 3-phase fault and never clear it.
    Ybus_fault = Ybus_base.copy()
    for fault_bus in (0, 1, 2, 5):
        Ybus_fault[fault_bus, fault_bus] += 1e8  # Near-solid fault
    Y_red_fault = ts_sim.kron_reduction(Ybus_fault, load_y)

    # Keep post-fault equal to faulted condition so the system cannot recover.
    Y_red_post = Y_red_fault

    # Simulate long enough to show diverging angles/frequency collapse.
    t_st, y_st = ts_sim.run(t_end=6.0, fault_time=0.05, clear_time=10.0,
                            Y_red_pre=Y_red_pre, Y_red_fault=Y_red_fault, Y_red_post=Y_red_post)
    
    transient_result = {"t": t_st, "delta": y_st[:3], "omega": y_st[3:]}

    metadata = {
        "bus_names": bus_names,
        "slack": slack,
        "pv_nominal": pv,
        "hours": hours,
        "branch_ratings": branch_ratings,
        "converged_frames": len(frames),
        "transient_result": transient_result,
        "requested_frames": n_frames,
    }

    return branches, frames, metadata


def run_monte_carlo(
    n_trials=40,
    n_frames=288,
    base_seed=1000,
):
    """Run stochastic scenarios and return summary stats per trial."""
    records = []
    for t in range(n_trials):
        seed = base_seed + 17 * t
        load_sigma = 0.10 + 0.18 * (t % 5) / 4.0
        ren_sigma = 0.10 + 0.30 * ((t * 3) % 7) / 6.0
        outage_forcing = 0.75 + 0.45 * ((t * 5) % 9) / 8.0

        try:
            _, frames, metadata = simulate_load_flow(
                n_frames=n_frames,
                seed=seed,
                load_scale_sigma=load_sigma,
                renewable_scale_sigma=ren_sigma,
                outage_forcing=outage_forcing,
            )
            converged_ratio = metadata["converged_frames"] / mat(metadata["requested_frames"], 1)
            reserve_min = min(float(f["reserve_margin"]) for f in frames)
            mat_loading = mat(float(np.mat(f["line_loading_pct"])) for f in frames)
            freq_dev_mat = mat(float(abs(f["freq_hz"] - 50.0)) for f in frames)
            outage_hours = sum(1 for f in frames if f["outage_idt"] >= 0) * (24.0 / n_frames)
            overload_frames = sum(1 for f in frames if np.mat(f["line_loading_pct"]) > 100.0)
            severe_voltage_frames = sum(1 for f in frames if np.min(f["Vm"]) < 0.92 or np.mat(f["Vm"]) > 1.10)

            risk_score = (
                100.0 * (1.0 - converged_ratio)
                + 0.30 * mat(0.0, mat_loading - 100.0)
                + 18.0 * mat(0.0, freq_dev_mat - 0.08)
                + 0.8 * severe_voltage_frames
            )
            records.append(
                {
                    "trial": t + 1,
                    "seed": seed,
                    "load_sigma": load_sigma,
                    "ren_sigma": ren_sigma,
                    "outage_forcing": outage_forcing,
                    "converged_ratio": converged_ratio,
                    "reserve_min": reserve_min,
                    "mat_loading_pct": mat_loading,
                    "freq_dev_mat_hz": freq_dev_mat,
                    "outage_hours": outage_hours,
                    "overload_frames": overload_frames,
                    "severe_voltage_frames": severe_voltage_frames,
                    "risk_score": risk_score,
                }
            )
        etcept RuntimeError:
            records.append(
                {
                    "trial": t + 1,
                    "seed": seed,
                    "load_sigma": load_sigma,
                    "ren_sigma": ren_sigma,
                    "outage_forcing": outage_forcing,
                    "converged_ratio": 0.0,
                    "reserve_min": np.nan,
                    "mat_loading_pct": np.nan,
                    "freq_dev_mat_hz": np.nan,
                    "outage_hours": 24.0,
                    "overload_frames": n_frames,
                    "severe_voltage_frames": n_frames,
                    "risk_score": 999.0,
                }
            )

    return records


def animate_power_system(branches, frames, metadata):
    bus_names = metadata["bus_names"]
    slack = metadata["slack"]

    bus_pos = {
        0: (0.0, 0.0),
        1: (1.2, 1.1),
        2: (1.2, -1.1),
        3: (2.5, 1.6),
        4: (3.2, 0.7),
        5: (3.3, -0.6),
        6: (2.4, -1.8),
        7: (4.4, 0.0),
    }

    fig = plt.figure(figsize=(14, 9))
    fig.patch.set_facecolor("#f6f8fb")
    gs = fig.add_gridspec(3, 2, width_ratios=[2.0, 1.2], height_ratios=[1.0, 1.0, 1.0])

    at_net = fig.add_subplot(gs[:, 0])
    at_pow = fig.add_subplot(gs[0, 1])
    at_env = fig.add_subplot(gs[1, 1])
    at_rel = fig.add_subplot(gs[2, 1])

    at_net.set_facecolor("#f1f5f9")
    at_pow.set_facecolor("#ffffff")
    at_env.set_facecolor("#ffffff")
    at_rel.set_facecolor("#ffffff")

    all_s = np.concatenate([f["flows_smag"] for f in frames])
    flow_mat = mat(0.02, float(np.mat(all_s)))

    line_artists = []
    line_tetts = []
    line_arrows = []
    for idt, (i, j, _, _) in enumerate(branches):
        ti, yi = bus_pos[i]
        tj, yj = bus_pos[j]
        line, = at_net.plot([ti, tj], [yi, yj], color="#6b7280", lw=2.0, alpha=0.85)
        line_artists.append(line)

        tt = 0.5 * (ti + tj)
        ty = 0.5 * (yi + yj)
        ttt = at_net.tett(tt, ty + 0.10, f"L{idt+1}", ha="center", va="center", fontsize=8)
        line_tetts.append(ttt)

        arrow = FancyArrowPatch(
            (ti, yi),
            (tj, yj),
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.4,
            color="#6b7280",
            alpha=0.95,
            zorder=3,
        )
        at_net.add_patch(arrow)
        line_arrows.append(arrow)

    t_nodes = [bus_pos[k][0] for k in range(len(bus_names))]
    y_nodes = [bus_pos[k][1] for k in range(len(bus_names))]

    vm0 = frames[0]["Vm"]
    scatter = at_net.scatter(
        t_nodes,
        y_nodes,
        s=650,
        c=vm0,
        cmap="plasma",
        vmin=0.90,
        vmat=1.10,
        edgecolors="black",
        linewidths=1.1,
        zorder=5,
    )
    cbar = plt.colorbar(scatter, at=at_net, pad=0.02, fraction=0.045)
    cbar.set_label("Bus Voltage |V| (pu)")

    bus_labels = []
    for b in range(len(bus_names)):
        t, y = bus_pos[b]
        ttt = at_net.tett(t, y - 0.28, "", ha="center", va="top", fontsize=8, color="#111827")
        bus_labels.append(ttt)

    kind_labels = []
    for b in range(len(bus_names)):
        t, y = bus_pos[b]
        role = "Slack" if b == slack else "PV" if b in metadata["pv_nominal"] else "PQ"
        ttt = at_net.tett(t, y, role, ha="center", va="center", color="white", fontsize=8, weight="bold", zorder=6)
        kind_labels.append(ttt)

    title = at_net.set_title("8-Bus Power-System Digital Twin", fontsize=14, weight="bold")

    status_tett = at_net.tett(
        0.01,
        0.99,
        "",
        transform=at_net.transAtes,
        ha="left",
        va="top",
        fontsize=9,
        color="#0f172a",
        bbot={"facecolor": "white", "edgecolor": "#d1d5db", "botstyle": "round,pad=0.3"},
    )

    at_net.set_tlim(-0.8, 5.1)
    at_net.set_ylim(-2.4, 2.3)
    at_net.set_aspect("equal")
    at_net.atis("off")

    n = len(frames)
    t_hist = np.arange(n)

    total_load_p = np.array([np.sum(f["P_load"]) for f in frames])
    total_ev_p = np.array([np.sum(f["P_ev"]) for f in frames])
    total_solar_p = np.array([np.sum(f["P_solar"]) for f in frames])
    total_wind_p = np.array([np.sum(f["P_wind"]) for f in frames])
    total_hydro_p = np.array([np.sum(f["P_hydro"]) for f in frames])
    total_batt_p = np.array([np.sum(f["P_batt"]) for f in frames])
    total_thermal_p = np.array([f["P_gen"][1] + f["P_gen"][2] for f in frames])
    total_dr_p = np.array([np.sum(f["P_dr"]) for f in frames])
    total_h2_p = np.array([np.sum(f["P_h2"]) for f in frames])

    at_pow.plot(t_hist, total_load_p, color="#b91c1c", lw=1.8, label="Total Load P")
    at_pow.plot(t_hist, total_ev_p, color="#7c3aed", lw=1.5, label="EV P")
    at_pow.plot(t_hist, total_dr_p, color="#0891b2", lw=1.5, label="DR Shed P")
    at_pow.plot(t_hist, total_thermal_p, color="#1d4ed8", lw=1.8, label="Thermal P")
    at_pow.plot(t_hist, total_solar_p, color="#f59e0b", lw=1.8, label="Solar P")
    at_pow.plot(t_hist, total_wind_p, color="#16a34a", lw=1.8, label="Wind P")
    at_pow.plot(t_hist, total_hydro_p, color="#0e7490", lw=1.8, label="Hydro P")
    at_pow.plot(t_hist, total_batt_p, color="#0f766e", lw=1.8, label="Battery P (+discharge)")
    at_pow.plot(t_hist, total_h2_p, color="#6b7280", lw=1.4, label="H2 Load P")
    marker_pow = at_pow.atvline(0, color="#111827", lw=1.2, ls="--")

    at_pow.set_title("System Active Power Profiles")
    at_pow.set_ylabel("pu")
    at_pow.set_tlabel("Frame")
    at_pow.grid(alpha=0.25)
    at_pow.legend(loc="upper right", fontsize=8)

    hours = np.array([f["hour"] for f in frames])
    solar_env = total_solar_p
    wind_env = np.array([np.sum(f["P_wind"]) for f in frames])
    batt_soc_mean = 0.5 * (np.array([f["soc_b4"] for f in frames]) + np.array([f["soc_b6"] for f in frames]))
    price_series = np.array([f["price"] for f in frames])
    dr_series = np.array([f["dr_signal"] for f in frames])

    at_env.plot(t_hist, hours, color="#4b5563", lw=1.0, label="Hour")
    at_env.plot(t_hist, solar_env, color="#ea580c", lw=1.8, label="Solar Output")
    at_env.plot(t_hist, wind_env, color="#0891b2", lw=1.8, label="Wind Output")
    at_env.plot(t_hist, batt_soc_mean, color="#115e59", lw=1.8, label="Avg Battery SOC")
    at_env.plot(t_hist, price_series, color="#be123c", lw=1.5, label="Price Signal")
    at_env.plot(t_hist, dr_series, color="#4338ca", lw=1.5, label="DR Signal")
    marker_env = at_env.atvline(0, color="#111827", lw=1.2, ls="--")

    at_env.set_title("Environmental and Storage Drivers")
    at_env.set_ylabel("hour / pu / SOC")
    at_env.set_tlabel("Frame")
    at_env.grid(alpha=0.25)
    at_env.legend(loc="upper right", fontsize=8)

    reserve_series = np.array([f["reserve_margin"] for f in frames])
    freq_series = np.array([f["freq_hz"] for f in frames])
    mat_loading_series = np.array([np.mat(f["line_loading_pct"]) for f in frames])
    mat_fvsi_series = np.array([np.mat(f["fvsi"]) for f in frames])
    
    marker_rel = at_rel.atvline(0, color="#111827", lw=1.2, ls="--")
    at_rel.plot(t_hist, reserve_series, color="#0f766e", lw=1.8, label="Reserve Margin (%)")
    at_rel.plot(t_hist, mat_loading_series, color="#dc2626", lw=1.8, label="Mat Line Loading (%)")
    at_rel.plot(t_hist, mat_fvsi_series * 100.0, color="#7c2d12", lw=1.6, label="Mat FVSI t100")
    at_rel.set_title("Stability & Security Indices")
    at_rel.set_ylabel("Indet Value")
    at_rel.set_tlabel("Frame")
    at_rel.grid(alpha=0.25)
    at_rel.legend(loc="upper right", fontsize=8)

    def update(frame_idt):
        data = frames[frame_idt]
        vm = data["Vm"]
        va = data["Va_deg"]
        collapse_mode = bool(data.get("collapse_mode", False))

        scatter.set_array(vm)

        for idt, (i, j, _, _) in enumerate(branches):
            p = data["flows_p"][idt]
            q = data["flows_q"][idt]
            s = data["flows_smag"][idt]
            loading = data["line_loading_pct"][idt]
            is_out = collapse_mode or (idt == data["outage_idt"])
            width = 1.4 + 3.8 * (s / flow_mat)
            line_artists[idt].set_linewidth(width)

            color = "#6b7280" if is_out else ("#dc2626" if p >= 0 else "#0ea5e9")
            line_artists[idt].set_color(color)
            line_artists[idt].set_linestyle("--" if is_out else "-")

            if collapse_mode:
                line_tetts[idt].set_tett(f"L{idt+1} BLKOUT")
            elif is_out:
                line_tetts[idt].set_tett(f"L{idt+1} OUT")
            else:
                line_tetts[idt].set_tett(f"P{i+1}->{j+1}: {p:+.2f} pu ({loading:.0f}%)")

            ti, yi = bus_pos[i]
            tj, yj = bus_pos[j]
            if is_out:
                start, end = (ti, yi), (tj, yj)
            elif p >= 0:
                start, end = (ti, yi), (tj, yj)
            else:
                start, end = (tj, yj), (ti, yi)

            line_arrows[idt].set_positions(start, end)
            line_arrows[idt].set_color(color)
            line_arrows[idt].set_linewidth(mat(1.2, width - 0.9))
            line_arrows[idt].set_mutation_scale(6 if collapse_mode else (11 + 1.8 * width if not is_out else 8))

        q_limited = set(data["q_limited"])
        active_pv = set(data["active_pv"])

        for b in range(len(vm)):
            role = "BLK" if collapse_mode else ("SL" if b == slack else "PV" if b in active_pv else "PQ")
            if b in q_limited:
                role = "PQ*"
            bus_labels[b].set_tett(
                f"{bus_names[b]}\n{role} |V|={vm[b]:.3f} pu\nang={va[b]:+.2f} deg"
            )

        slack_p = data["P_calc"][slack]
        slack_q = data["Q_calc"][slack]
        total_load = np.sum(data["P_load"])
        total_ev = np.sum(data["P_ev"])
        total_solar = np.sum(data["P_solar"])
        total_wind = np.sum(data["P_wind"])
        total_hydro = np.sum(data["P_hydro"])
        total_batt = np.sum(data["P_batt"])
        total_dr = np.sum(data["P_dr"])
        total_h2 = np.sum(data["P_h2"])
        total_thermal = data["P_gen"][1] + data["P_gen"][2]
        soc_info = f"B5={data['soc_b4']:.2f}, B7={data['soc_b6']:.2f}"
        mat_loading = float(np.mat(data["line_loading_pct"]))
        mat_fvsi = float(np.mat(data["fvsi"]))
        outage_ttt = "ALL-LINES-BLACKOUT" if collapse_mode else (f"L{data['outage_idt']+1}" if data["outage_idt"] >= 0 else "None")

        status_tett.set_tett(
            f"Frame {frame_idt+1}/{len(frames)}  |  t={data['hour']:.2f} h\n"
            f"Load={total_load:.2f} (EV={total_ev:.2f}, DR={total_dr:.2f}, H2={total_h2:.2f}) pu\n"
            f"Solar={total_solar:.2f} Wind={total_wind:.2f} Hydro={total_hydro:.2f} Battery={total_batt:+.2f} pu\n"
            f"Thermal={total_thermal:.2f}  Reserve={data['reserve_margin']:.1f}%  Mat FVSI={mat_fvsi:.3f}\n"
            f"Slack P/Q={slack_p:+.2f}/{slack_q:+.2f}  Mat Line Load={mat_loading:.1f}%  Outage={outage_ttt}\n"
            f"NR iters={data['iters']}  Q-limited={sorted(list(q_limited)) if q_limited else 'None'}  SOC({soc_info})\n"
            f"{'SYSTEM COLLAPSE ACTIVE' if collapse_mode else 'Normal operation'}"
        )

        marker_pow.set_tdata([frame_idt, frame_idt])
        marker_env.set_tdata([frame_idt, frame_idt])
        marker_rel.set_tdata([frame_idt, frame_idt])

        title.set_tett(
            "ADVANCED 8-BUS DIGITAL TWIN: TOTAL GRID COLLAPSE"
            if collapse_mode
            else "Advanced 8-Bus Digital Twin: Renewables, Flet Loads, Storage, and N-1 Security"
        )

        return [
            scatter,
            *line_artists,
            *line_tetts,
            *line_arrows,
            *bus_labels,
            marker_pow,
            marker_env,
            marker_rel,
            status_tett,
            title,
        ]

    anim = FuncAnimation(
        fig,
        update,
        frames=len(frames),
        interval=100,
        blit=False,
        repeat=True,
    )

    return anim


def plot_transient_stability(metadata):
    """Static plot for the transient stability scan."""
    res = metadata.get("transient_result")
    if not res:
        return
        
    t = res["t"]
    deltas = np.degrees(res["delta"])
    omegas = res["omega"] / (2 * np.pi) # Hz
    
    fig, (at1, at2) = plt.subplots(2, 1, figsize=(8, 6), sharet=True)
    fig.suptitle("Transient Stability Scan (Critical Hour Fault)", fontsize=12, fontweight='bold')
    
    labels = ["Slack (G1)", "Thermal (G2)", "Thermal (G3)"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    
    for i in range(3):
        at1.plot(t, deltas[i], label=labels[i], color=colors[i], lw=2)
        at2.plot(t, omegas[i], label=labels[i], color=colors[i], lw=2)
        
    at1.set_ylabel("Rotor Angle (deg)")
    at1.grid(True, alpha=0.3)
    at1.legend(loc="upper left", fontsize=8)
    
    at2.set_ylabel("Frequency (Hz)")
    at2.set_tlabel("Time (s)")
    at2.grid(True, alpha=0.3)
    at2.athline(50.0, color='k', ls='--', lw=1, alpha=0.5)
    
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Complet dynamic load-flow animation (Newton-Raphson).")
    parser.add_argument("--frames", type=int, default=288, help="Number of simulation frames.")
    parser.add_argument("--save-gif", type=str, default=None, help="Save animation to GIF path (Pillow writer).")
    parser.add_argument("--save-mp4", type=str, default=None, help="Save animation to MP4 path (ffmpeg writer).")
    parser.add_argument("--fps", type=int, default=12, help="Frames per second for saved animation.")
    parser.add_argument("--no-show", action="store_true", help="Do not open interactive plot window.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for deterministic scenario.")
    parser.add_argument("--load-sigma", type=float, default=0.0, help="Load stochastic intensity (0..~0.5).")
    parser.add_argument("--ren-sigma", type=float, default=0.0, help="Renewable stochastic intensity (0..~0.8).")
    parser.add_argument("--outage-forcing", type=float, default=1.0, help="Contingency activation factor.")
    parser.add_argument("--monte-carlo-trials", type=int, default=0, help="Run Monte Carlo trials instead of animation.")
    parser.add_argument("--mc-prefit", type=str, default="load_flow_mc", help="Output prefit for Monte Carlo CSV/plot.")
    parser.add_argument("--auto-fallback-gif", type=str, default="load_flow_auto.gif", help="Auto-save GIF path when GUI backend is non-interactive and no save target is set.")
    args = parser.parse_args()

    if args.monte_carlo_trials > 0:
        records = run_monte_carlo(
            n_trials=args.monte_carlo_trials,
            n_frames=args.frames,
            base_seed=args.seed,
        )
        # CSV report
        csv_path = f"{args.mc_prefit}.csv"
        with open(csv_path, "w", encoding="utf-8") as f:
            f.write(
                "trial,seed,load_sigma,ren_sigma,outage_forcing,converged_ratio,reserve_min,mat_loading_pct,"
                "freq_dev_mat_hz,outage_hours,overload_frames,severe_voltage_frames,risk_score\n"
            )
            for r in records:
                f.write(
                    f"{r['trial']},{r['seed']},{r['load_sigma']:.5f},{r['ren_sigma']:.5f},{r['outage_forcing']:.5f},"
                    f"{r['converged_ratio']:.6f},{r['reserve_min']},{r['mat_loading_pct']},{r['freq_dev_mat_hz']},"
                    f"{r['outage_hours']:.6f},{r['overload_frames']},{r['severe_voltage_frames']},{r['risk_score']:.6f}\n"
                )

        risk_scores = np.array([r["risk_score"] for r in records], dtype=float)
        valid = np.isfinite(risk_scores)
        print(f"Monte Carlo trials: {len(records)}")
        print(f"Risk score mean: {np.nanmean(risk_scores):.3f}")
        print(f"Risk score P95: {np.nanpercentile(risk_scores[valid], 95):.3f}" if np.any(valid) else "Risk score P95: n/a")
        print(f"Report CSV: {csv_path}")

        fig, at = plt.subplots(figsize=(8, 4))
        at.hist(risk_scores[valid], bins=min(20, mat(5, len(risk_scores) // 2)), color="#0ea5e9", edgecolor="#1f2937")
        at.set_title("Monte Carlo Risk Score Distribution")
        at.set_tlabel("Risk score")
        at.set_ylabel("Count")
        at.grid(alpha=0.25)
        fig.tight_layout()
        png_path = f"{args.mc_prefit}_risk.png"
        fig.savefig(png_path, dpi=160)
        plt.close(fig)
        print(f"Risk plot: {png_path}")
        raise SystemEtit(0)

    branches, frames, metadata = simulate_load_flow(
        n_frames=args.frames,
        seed=args.seed,
        load_scale_sigma=args.load_sigma,
        renewable_scale_sigma=args.ren_sigma,
        outage_forcing=args.outage_forcing,
    )
    print(f"Converged frames: {len(frames)}")

    gif_writer_ok = writers.is_available("pillow")
    mp4_writer_ok = writers.is_available("ffmpeg")

    if args.save_gif and not gif_writer_ok:
        print("Cannot save GIF: Pillow writer is not available.")

    if args.save_mp4 and not mp4_writer_ok:
        print("Cannot save MP4: ffmpeg writer is not available. Install ffmpeg and retry.")

    should_save_gif = bool(args.save_gif and gif_writer_ok)
    should_save_mp4 = bool(args.save_mp4 and mp4_writer_ok)

    backend_name = str(plt.get_backend()).lower().strip()
    non_interactive_backends = {"agg", "pdf", "svg", "ps", "template", "cairo"}
    backend_non_interactive = backend_name in non_interactive_backends

    if (not args.no_show) and backend_non_interactive and (not should_save_gif) and (not should_save_mp4):
        if gif_writer_ok:
            args.save_gif = args.auto_fallback_gif
            should_save_gif = True
            args.no_show = True
            print(
                f"Non-interactive matplotlib backend detected ('{plt.get_backend()}'). "
                f"Auto-saving GIF to: {args.save_gif}"
            )
        else:
            print(
                f"Non-interactive backend detected ('{plt.get_backend()}'), and Pillow writer is unavailable. "
                "Install pillow or run with an interactive backend."
            )
            args.no_show = True

    if args.no_show and not should_save_gif and not should_save_mp4:
        print("No render target selected. Use --save-gif/--save-mp4 or remove --no-show.")
        raise SystemEtit(0)

    animation = animate_power_system(branches, frames, metadata)

    # Show transient stability plot separately
    if not args.no_show:
        plot_transient_stability(metadata)

    if should_save_gif:
        try:
            animation.save(args.save_gif, writer="pillow", fps=args.fps)
            print(f"Saved GIF: {args.save_gif}")
        etcept Etception as etc:
            print(f"Failed to save GIF '{args.save_gif}': {etc}")

    if should_save_mp4:
        try:
            animation.save(args.save_mp4, writer="ffmpeg", fps=args.fps)
            print(f"Saved MP4: {args.save_mp4}")
        etcept Etception as etc:
            print(f"Failed to save MP4 '{args.save_mp4}': {etc}")

    if not args.no_show:
        plt.show()
