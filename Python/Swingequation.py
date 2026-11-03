"""
Power System Transient Stability Simulation — RK4 Method
==========================================================
Model: Single Machine Infinite Bus (SMIB) with classical machine model.

Swing Equation:
    M · d²δ/dt² = Pm − Pe − D·(dδ/dt)

Rewritten as a first-order system:
    dδ/dt   = Δω
    dΔω/dt  = (π·f / H) · (Pm − Pe − D·Δω)

Electrical Power:
    Pe = (E' · V∞ / X_total) · sin(δ)

Three topological states are simulated:
    • Pre-fault   — both lines in service
    • During fault— 3-phase bolted fault at midpoint of one line
    • Post-fault  — faulty line tripped, single line remaining

Critical clearing time is estimated via the Equal Area Criterion (EAC).
"""

import numpy as np
import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.collections import LineCollection
from scipy.integrate import trapezoid



# ──────────────────────────────────────────────────────────────────
class PowerSystemSimulator:
    # ── constructor ──────────────────────────────────────────────
    def __init__(self, root):
        self.root = root
        self.root.title("Power System Transient Stability — RK4 Simulation")
        self.root.configure(bg="#f0f0f0")

        # ---------- default electrical / mechanical parameters ----------
        self.H          = tk.DoubleVar(value=5.0)    # inertia constant  [s]
        self.D          = tk.DoubleVar(value=2.0)    # damping coeff.    [pu]
        self.Xd_prime   = tk.DoubleVar(value=0.30)   # d-axis transient reactance [pu]
        self.Xt         = tk.DoubleVar(value=0.10)   # transformer reactance     [pu]
        self.X_line     = tk.DoubleVar(value=0.40)   # line reactance per circuit [pu]
        self.E_prime    = tk.DoubleVar(value=1.20)   # internal voltage behind Xd' [pu]
        self.V_inf      = tk.DoubleVar(value=1.00)   # infinite-bus voltage [pu]
        self.Pm         = tk.DoubleVar(value=0.80)   # mechanical power    [pu]

        # ---------- fault / simulation timing ----------
        self.fault_time = tk.DoubleVar(value=0.10)   # fault inception     [s]
        self.clear_time = tk.DoubleVar(value=0.30)   # fault clearing      [s]
        self.sim_time   = tk.DoubleVar(value=2.00)   # total sim duration  [s]

        # ---------- numerical settings ----------
        self.dt  = 0.0005          # RK4 time-step [s]
        self.f0  = 60.0            # nominal frequency [Hz]
        self.cbar = None           # colorbar handle

        # value-label lookup (populated in build_controls)
        self._val_labels = {}

        self._build_gui()
        self._run_simulation()

    # ── GUI layout ───────────────────────────────────────────────
    def _build_gui(self):
        style = ttk.Style()
        style.configure("Header.TLabel", font=("Segoe UI", 9, "bold"))

        outer = ttk.Frame(self.root)
        outer.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)

        # ---- left: controls ----
        ctrl = ttk.LabelFrame(outer, text="  Parameters  ", padding=8)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 6))

        slider_defs = [
            ("H  (s)",       self.H,          1.0,  10.0,  0.5),
            ("D  (pu)",      self.D,          0.0,  10.0,  0.5),
            ("Xd' (pu)",    self.Xd_prime,   0.10,  0.50, 0.02),
            ("Xt (pu)",     self.Xt,         0.00,  0.30, 0.01),
            ("Xline (pu)",  self.X_line,     0.10,  1.00, 0.05),
            ("E' (pu)",     self.E_prime,    0.80,  1.50, 0.05),
            ("V∞ (pu)",     self.V_inf,      0.80,  1.20, 0.05),
            ("Pm (pu)",     self.Pm,         0.10,  1.50, 0.05),
            ("tf (s)",      self.fault_time, 0.00,  1.00, 0.01),
            ("tc (s)",      self.clear_time, 0.02,  1.00, 0.01),
            ("Tsim (s)",    self.sim_time,   0.5,   5.0,  0.1),
        ]

        for row, (label, var, lo, hi, res) in enumerate(slider_defs):
            ttk.Label(ctrl, text=label, style="Header.TLabel").grid(
                row=row, column=0, sticky="w", pady=3)
            sc = tk.Scale(ctrl, variable=var, from_=lo, to=hi,
                          resolution=res, orient=tk.HORIZONTAL, length=140,
                          showvalue=False, command=self._on_slider)
            sc.grid(row=row, column=1, padx=4, pady=2)
            vl = ttk.Label(ctrl, text=f"{var.get():.3f}", width=7, anchor="e")
            vl.grid(row=row, column=2, pady=2)
            self._val_labels[label] = vl

        sep = ttk.Separator(ctrl, orient=tk.HORIZONTAL)
        sep.grid(row=len(slider_defs), column=0, columnspan=3,
                 sticky="ew", pady=8)

        btn_row = len(slider_defs) + 1
        ttk.Button(ctrl, text="▶  Run Simulation",
                   command=self._run_simulation).grid(
            row=btn_row, column=0, columnspan=3, pady=4, sticky="ew")
        ttk.Button(ctrl, text="↺  Reset Defaults",
                   command=self._reset_defaults).grid(
            row=btn_row+1, column=0, columnspan=3, pady=4, sticky="ew")

        self.status_var = tk.StringVar()
        ttk.Label(ctrl, textvariable=self.status_var, wraplength=210,
                  foreground="#0055aa", font=("Segoe UI", 9, "bold")).grid(
            row=btn_row+2, column=0, columnspan=3, pady=(10, 2), sticky="w")

        self.info_var = tk.StringVar()
        ttk.Label(ctrl, textvariable=self.info_var, wraplength=220,
                  justify="left", font=("Consolas", 8)).grid(
            row=btn_row+3, column=0, columnspan=3, pady=4, sticky="w")

        # ---- right: matplotlib canvas ----
        plot_frame = ttk.Frame(outer)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.fig, self.axes = plt.subplots(2, 2, figsize=(12, 8),
                                           facecolor="#fafafa")
        self.fig.subplots_adjust(left=0.07, right=0.95, bottom=0.08,
                                 top=0.94, wspace=0.34, hspace=0.38)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # ── helpers ──────────────────────────────────────────────────
    def _on_slider(self, *_):
        for label, var in [("H  (s)", self.H), ("D  (pu)", self.D),
                           ("Xd' (pu)", self.Xd_prime), ("Xt (pu)", self.Xt),
                           ("Xline (pu)", self.X_line), ("E' (pu)", self.E_prime),
                           ("V∞ (pu)", self.V_inf), ("Pm (pu)", self.Pm),
                           ("tf (s)", self.fault_time), ("tc (s)", self.clear_time),
                           ("Tsim (s)", self.sim_time)]:
            self._val_labels[label].config(text=f"{var.get():.3f}")

    def _reset_defaults(self):
        defaults = {self.H: 5.0, self.D: 2.0, self.Xd_prime: 0.30,
                    self.Xt: 0.10, self.X_line: 0.40, self.E_prime: 1.20,
                    self.V_inf: 1.00, self.Pm: 0.80,
                    self.fault_time: 0.10, self.clear_time: 0.30,
                    self.sim_time: 2.00}
        for var, val in defaults.items():
            var.set(val)
        self._on_slider()
        self._run_simulation()

    # ── electrical model ─────────────────────────────────────────
    def _reactances(self):
        """Return (X_pre, X_fault, X_post) in per-unit."""
        Xd = self.Xd_prime.get()
        Xt = self.Xt.get()
        Xl = self.X_line.get()
        X_pre   = Xd + Xt + Xl / 2.0          # two lines in parallel
        X_fault = Xd + Xt + Xl + Xl            # fault at midpoint → very high Z
        X_post  = Xd + Xt + Xl                 # one line tripped
        return X_pre, X_fault, X_post

    def _Pe(self, delta, X):
        return self.E_prime.get() * self.V_inf.get() / X * np.sin(delta)

    # ── differential equations ───────────────────────────────────
    def _derivs(self, state, t):
        """
        state = [δ, Δω]   (Δω = dδ/dt  in rad/s)

        dδ/dt   = Δω
        dΔω/dt  = (π f₀ / H) · (Pm − Pe − D Δω)
        """
        delta, omega = state
        X_pre, X_fault, X_post = self._reactances()
        tf, tc = self.fault_time.get(), self.clear_time.get()

        X = X_pre if t < tf else (X_fault if t < tc else X_post)
        Pe = self._Pe(delta, X)

        ddelta = omega
        domega = (np.pi * self.f0 / self.H.get()) * \
                 (self.Pm.get() - Pe - self.D.get() * omega)
        return np.array([ddelta, domega])

    # ── RK4 integrator ──────────────────────────────────────────
    def _rk4_step(self, y, t, h):
        k1 = self._derivs(y, t)
        k2 = self._derivs(y + 0.5 * h * k1, t + 0.5 * h)
        k3 = self._derivs(y + 0.5 * h * k2, t + 0.5 * h)
        k4 = self._derivs(y + h * k3, t + h)
        return y + (h / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)

    # ── critical clearing time via EAC ──────────────────────────
    def _critical_clearing_time(self, delta0, delta_crit):
        E, V = self.E_prime.get(), self.V_inf.get()
        X_pre, X_fault, X_post = self._reactances()
        Pm = self.Pm.get()
        Pmax_f = E * V / X_fault
        Pmax_p = E * V / X_post
        H, f0 = self.H.get(), self.f0

        # --- find δ_clear where A1 = A2 ---
        best_d, best_err = delta0, 1e9
        dc = delta0 + 1e-4
        while dc < delta_crit - 1e-4:
            d1 = np.linspace(delta0, dc, 200)
            A1 = trapezoid(Pm - Pmax_f * np.sin(d1), d1)
            d2 = np.linspace(dc, delta_crit, 200)
            A2 = trapezoid(Pmax_p * np.sin(d2) - Pm, d2)
            err = abs(A1 - A2)
            if err < best_err:
                best_err, best_d = err, dc
            dc += 1e-4

        # --- integrate swing eq. during fault to find time to best_d ---
        dt_fine = 1e-5
        state = np.array([delta0, 0.0])
        t = 0.0
        while state[0] < best_d and t < 5.0:
            # fault-on derivatives (no damping for EAC)
            d, w = state
            def fault_derivs(s, _t):
                dd = s[1]
                dw = (np.pi * f0 / H) * (Pm - Pmax_f * np.sin(s[0]))
                return np.array([dd, dw])
            k1 = fault_derivs(state, t)
            k2 = fault_derivs(state + 0.5*dt_fine*k1, t+0.5*dt_fine)
            k3 = fault_derivs(state + 0.5*dt_fine*k2, t+0.5*dt_fine)
            k4 = fault_derivs(state + dt_fine*k3, t+dt_fine)
            state = state + (dt_fine/6.0)*(k1 + 2*k2 + 2*k3 + k4)
            t += dt_fine
        return t, best_d

    # ── main simulation loop ────────────────────────────────────
    def _run_simulation(self):
        self._on_slider()

        X_pre, X_fault, X_post = self._reactances()
        E, V = self.E_prime.get(), self.V_inf.get()
        Pm   = self.Pm.get()
        tf   = self.fault_time.get()
        tc   = self.clear_time.get()
        T    = self.sim_time.get()
        dt   = self.dt

        Pmax_pre = E * V / X_pre
        if Pm > Pmax_pre:
            self.status_var.set("⚠ No pre-fault equilibrium (Pm > Pmax)!")
            return

        # initial steady-state
        delta0 = np.arcsin(Pm / Pmax_pre)
        omega0 = 0.0

        # allocate
        N = int(T / dt) + 1
        t_arr     = np.linspace(0, T, N)
        delta_arr = np.zeros(N)
        omega_arr = np.zeros(N)
        Pe_arr    = np.zeros(N)

        state = np.array([delta0, omega0])
        delta_arr[0], omega_arr[0] = state
        Pe_arr[0] = self._Pe(delta0, X_pre)

        # ── RK4 integration ──
        for i in range(1, N):
            state = self._rk4_step(state, t_arr[i-1], dt)
            delta_arr[i] = state[0]
            omega_arr[i] = state[1]
            ti = t_arr[i]
            X = X_pre if ti < tf else (X_fault if ti < tc else X_post)
            Pe_arr[i] = self._Pe(state[0], X)

        # derived arrays
        delta_deg = np.degrees(delta_arr)
        omega_pu  = omega_arr / (2.0 * np.pi * self.f0)   # per-unit speed dev

        # stability verdict (look at last 0.5 s)
        tail = max(1, int(0.5 / dt))
        stable = (np.max(np.abs(delta_deg[-tail:])) < 180 and
                  np.max(np.abs(omega_pu[-tail:])) < 0.05)

        # critical clearing info
        delta_crit = np.pi - np.arcsin(Pm / (E * V / X_post))
        tcc, dcc = self._critical_clearing_time(delta0, delta_crit)

        # actual clearing angle from simulation
        idx_c = min(int(tc / dt), N - 1)
        delta_clear_sim = delta_arr[idx_c]

        # ── info panel ──
        info = (
            f"── Reactances ──────────\n"
            f" X_pre    = {X_pre:.4f} pu\n"
            f" X_fault  = {X_fault:.4f} pu\n"
            f" X_post   = {X_post:.4f} pu\n"
            f"── Pmax values ─────────\n"
            f" Pre      = {E*V/X_pre:.3f} pu\n"
            f" Fault    = {E*V/X_fault:.3f} pu\n"
            f" Post     = {E*V/X_post:.3f} pu\n"
            f"── Key angles ──────────\n"
            f" δ₀       = {np.degrees(delta0):.2f}°\n"
            f" δ_clear  = {np.degrees(delta_clear_sim):.2f}°\n"
            f" δ_crit   = {np.degrees(delta_crit):.2f}°\n"
            f"── Clearing times ──────\n"
            f" t_crit   ≈ {tcc:.4f} s\n"
            f" t_clear  = {tc:.4f} s\n"
            f" margin   = {(tc-tcc)*1000:.1f} ms"
        )
        self.info_var.set(info)

        if stable:
            self.status_var.set(
                f"✓ STABLE\n"
                f"δ_final = {delta_deg[-1]:.2f}°  "
                f"Δω_final = {omega_pu[-1]*100:.4f}%")
        else:
            self.status_var.set(
                f"✗ UNSTABLE — loss of synchronism\n"
                f"δ_max = {np.max(delta_deg):.1f}°")

        # ── update plots ──
        self._draw(t_arr, delta_deg, omega_pu, Pe_arr,
                   Pm, X_pre, X_fault, X_post,
                   tf, tc, delta0, delta_clear_sim, delta_crit)

    # ── plotting ─────────────────────────────────────────────────
    def _draw(self, t, delta_deg, omega_pu, Pe, Pm,
              X_pre, X_fault, X_post,
              tf, tc, delta0, delta_clear, delta_crit):

        # remove old colorbar if present
        if self.cbar is not None:
            try:
                self.cbar.remove()
            except Exception:
                pass
            self.cbar = None

        E, V = self.E_prime.get(), self.V_inf.get()
        Pmax_pre   = E * V / X_pre
        Pmax_fault = E * V / X_fault
        Pmax_post  = E * V / X_post

        for ax in self.axes.flat:
            ax.clear()

        # ── 1  Swing curve ──────────────────────────────────────
        ax1 = self.axes[0, 0]
        ax1.plot(t, delta_deg, color="#1f77b4", lw=1.6)
        ax1.axvline(tf, color="red",   ls="--", lw=1, alpha=.7, label=f"Fault @ {tf:.3f} s")
        ax1.axvline(tc, color="green", ls="--", lw=1, alpha=.7, label=f"Clear @ {tc:.3f} s")
        ax1.axhline(180, color="k", ls=":", lw=.8, alpha=.4, label="180° limit")
        ax1.axhline(np.degrees(delta0), color="gray", ls="-.", lw=.8,
                    alpha=.5, label=f"δ₀ = {np.degrees(delta0):.1f}°")
        ax1.axvspan(tf, tc, color="red", alpha=0.08)
        ax1.set_xlabel("Time  (s)")
        ax1.set_ylabel("Rotor Angle  δ  (°)")
        ax1.set_title("Swing Curve", fontweight="bold")
        ax1.legend(fontsize=7, loc="best")
        ax1.grid(True, alpha=0.25)
        ax1.set_xlim(0, t[-1])

        # ── 2  Phase portrait ───────────────────────────────────
        ax2 = self.axes[0, 1]
        step = max(1, len(t) // 3000)
        ds, ws, ts = delta_deg[::step], omega_pu[::step], t[::step]
        pts  = np.column_stack([ds, ws]).reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        lc   = LineCollection(segs, cmap="plasma", lw=1.6)
        lc.set_array(ts[:-1])
        ax2.add_collection(lc)
        self.cbar = self.fig.colorbar(lc, ax=ax2, shrink=0.82,
                                       label="Time (s)", pad=0.02)

        # equilibria
        ax2.plot(np.degrees(delta0), 0, "go", ms=10, zorder=5, label="Pre-fault eq.")
        d_post_eq = np.arcsin(Pm / Pmax_post)
        ax2.plot(np.degrees(d_post_eq), 0, "rs", ms=10, zorder=5,
                 label="Post-fault eq. (stable)")
        ax2.plot(np.degrees(np.pi - d_post_eq), 0, "r^", ms=10, zorder=5,
                 label="Post-fault eq. (unstable)")
        ax2.plot(delta_deg[0], omega_pu[0], "k*", ms=14, zorder=6, label="Start")

        ax2.autoscale_view()
        ax2.set_xlabel("Rotor Angle  δ  (°)")
        ax2.set_ylabel("Speed Deviation  Δω  (pu)")
        ax2.set_title("Phase Portrait", fontweight="bold")
        ax2.legend(fontsize=6.5, loc="best")
        ax2.grid(True, alpha=0.25)

        # ── 3  Power-angle + EAC ───────────────────────────────
        ax3 = self.axes[1, 0]
        d = np.linspace(0, np.pi, 600)
        ax3.plot(np.degrees(d), Pmax_pre   * np.sin(d), "b-",  lw=2,
                 label=f"Pre  (Pmax={Pmax_pre:.2f})")
        ax3.plot(np.degrees(d), Pmax_fault * np.sin(d), "r-",  lw=2,
                 label=f"Fault (Pmax={Pmax_fault:.2f})")
        ax3.plot(np.degrees(d), Pmax_post  * np.sin(d), "g-",  lw=2,
                 label=f"Post (Pmax={Pmax_post:.2f})")
        ax3.axhline(Pm, color="k", ls="--", lw=1.3, label=f"Pm = {Pm:.2f}")

        # accelerating area A1
        d1 = np.linspace(delta0, delta_clear, 150)
        p1 = Pmax_fault * np.sin(d1)
        ax3.fill_between(np.degrees(d1), Pm, p1,
                         where=(Pm >= p1), color="red", alpha=0.30,
                         label="A₁ (accel.)")
        # decelerating area A2
        d2 = np.linspace(delta_clear, delta_crit, 150)
        p2 = Pmax_post * np.sin(d2)
        ax3.fill_between(np.degrees(d2), Pm, p2,
                         where=(p2 >= Pm), color="green", alpha=0.30,
                         label="A₂ (decel.)")

        for ang, c, lbl in [(delta0, "blue", "δ₀"),
                             (delta_clear, "red", "δ_cl"),
                             (delta_crit, "green", "δ_cr")]:
            ax3.axvline(np.degrees(ang), color=c, ls=":", lw=.9, alpha=.6)
            ax3.text(np.degrees(ang), Pmax_pre * 1.05, lbl,
                     ha="center", fontsize=7, color=c)

        ax3.set_xlim(0, 180)
        ax3.set_ylim(0, max(Pmax_pre, Pmax_post) * 1.15)
        ax3.set_xlabel("Rotor Angle  δ  (°)")
        ax3.set_ylabel("Power  (pu)")
        ax3.set_title("Power-Angle Curves & Equal Area Criterion",
                      fontweight="bold")
        ax3.legend(fontsize=6, loc="upper right")
        ax3.grid(True, alpha=0.25)

        # ── 4  Pe(t) and Δω(t) ─────────────────────────────────
        ax4 = self.axes[1, 1]
        c1 = "#1f77b4"
        ax4.plot(t, Pe, color=c1, lw=1.4, label="Pe(t)")
        ax4.axhline(Pm, color="k", ls="--", lw=1, alpha=.6,
                    label=f"Pm = {Pm:.2f}")
        ax4.axvline(tf, color="red",   ls="--", lw=.8, alpha=.5)
        ax4.axvline(tc, color="green", ls="--", lw=.8, alpha=.5)
        ax4.axvspan(tf, tc, color="red", alpha=0.08)
        ax4.set_xlabel("Time  (s)")
        ax4.set_ylabel("Electrical Power  Pe  (pu)", color=c1)
        ax4.tick_params(axis="y", labelcolor=c1)
        ax4.set_xlim(0, t[-1])

        ax4b = ax4.twinx()
        c2 = "#d62728"
        ax4b.plot(t, omega_pu * 100, color=c2, lw=1.0, alpha=.85,
                  label="Δω(t)")
        ax4b.set_ylabel("Speed Deviation  Δω  (%)", color=c2)
        ax4b.tick_params(axis="y", labelcolor=c2)

        h1, l1 = ax4.get_legend_handles_labels()
        h2, l2 = ax4b.get_legend_handles_labels()
        ax4.legend(h1 + h2, l1 + l2, fontsize=7, loc="best")
        ax4.set_title("Electrical Power & Speed Deviation",
                      fontweight="bold")
        ax4.grid(True, alpha=0.25)

        self.canvas.draw_idle()


# ══════════════════════════════════════════════════════════════════
def main():
    root = tk.Tk()
    root.geometry("1450x880")
    root.minsize(1100, 700)
    PowerSystemSimulator(root)
    root.mainloop()


if __name__ == "__main__":
    main()