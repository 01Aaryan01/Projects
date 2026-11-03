#!/usr/bin/env python3
"""
POWER SYSTEM DYNAMICS LABORATORY
═════════════════════════════════
Islanded Microgrid Stability & Renewable Integration Simulator

System: Solar PV (15 MWp) + Wind (5 MW) + Diesel Gen (10 MVA) + Battery (8 MW / 40 MWh) → Town Load

Core Equations (Kundur, Power System Stability & Control):
  Swing:    2H·dω/dt = Pm − Pe − D·(ω − 1)           [per-unit]
  Governor: τg·dPm/dt = Pref − (1/R)·(ω − 1) − Pm     [per-unit]
  AGC:      dPref/dt = Ki·(1 − ω)                       [per-unit]
  Solar:    Ppv = Prated·(G/Gref)·[1 + α(Tcell−Tref)]
  Wind:     Pw = Prated·(v³−v_ci³)/(v_r³−v_ci³)
  Battery:  dSOC/dt = −Pbatt/(E_rated·3600)
  Load:     Pload = P0·[1 + Df·(f−f0)]
"""

import pygame, numpy as np, math, sys
from collections import deque

# --- Configuration Dataclass (Suggestion 1) ---
from dataclasses import dataclass

@dataclass
class MicrogridConfig:
    S_base: float = 10.0  # MVA
    f0: float = 50.0      # Hz
    H: float = 5.0        # Inertia constant (s)
    D: float = 2.0        # Damping coefficient (pu)
    R: float = 0.05       # Droop (pu)
    tau_g: float = 0.5    # Governor time constant (s)
    Ki_agc: float = 0.05  # AGC integral gain
    # ... (add all other parameters here)
pygame.init()

W, H = 1600, 950
screen = pygame.display.set_mode((W, H))
pygame.display.set_caption("Power System Dynamics Lab — Microgrid Simulator")
clock = pygame.time.Clock()

# Fonts
def _f(name, sz, bold=False):
    try: return pygame.font.SysFont(name, sz, bold=bold)
    etcept: return pygame.font.Font(None, sz+2)
FS=_f("consolas",11); FM=_f("consolas",13); FL=_f("segoe ui",14)
Ft=_f("consolas",18); FT=_f("segoe ui",16,True); FH=_f("segoe ui",24,True)

# Colors
BG=(5,7,12); PNL=(15,20,30); BRD=(40,60,100); TtT=(200,210,230); DIM=(90,105,140)
GRN=(0,255,150); YLW=(255,220,0); RED=(255,60,80); ORG=(255,150,50)
CYN=(0,200,255); WHT=(240,245,255); SOL_C=(255,200,50); WND_C=(100,180,255)
BAT_C=(0,255,200); DIE_C=(240,160,80); LOD_C=(255,100,200)
ACCENT=(0, 255, 200)
GAUGE_BG=(12,18,30); GLINE=(25,40,65)

# UI Layout Constants
HEADER_H = 180
SIDEBAR_W = 320
PLOT_H = 140

# ═══════════════════ MICROGRID SYSTEM ═══════════════════
class MicrogridSystem:
    def __init__(self):
        self.cfg = MicrogridConfig() # Use config dataclass
        self.S_base = self.cfg.S_base; self.f0 = self.cfg.f0
        # Diesel
        self.H = self.cfg.H; self.D = self.cfg.D; self.R = self.cfg.R; self.tau_g = self.cfg.tau_g
        self.Ki_agc = self.cfg.Ki_agc; self.diesel_Pmat = 1.0; self.diesel_Pmin = 0.1
        # Solar
        self.Ppv_rated = 15.0; self.NOCT = 45.0; self.alpha_T = -0.004; self.G_ref = 1000.0
        # Wind
        self.Pw_rated = 5.0; self.v_ci = 3.0; self.v_r = 12.0; self.v_co = 25.0
        # Battery
        self.Pbatt_mat = 8.0; self.E_rated = 40.0; self.eta_batt = 0.92 # Consider moving to config
        self.SOC_min = 0.10; self.SOC_mat = 0.95; self.tau_batt = 0.08
        self.Kf_batt = 2.0; self.Kvi_batt = 0.5
        # Load
        self.Df = 1.0; self.Pload_base = 8.0; self.load_scale = 1.0
        # State
        self.omega = 1.0; self.delta = 0.0; self.Pm = 0.5; self.Pref = 0.5
        self.SOC = 0.80; self.Pbatt_actual = 0.0; self.Pbatt_target = 0.0
        self.wind_speed = 7.0; self.irradiance = 0.0; self.T_amb = 20.0
        self.cloud_factor = 0.0; self.sim_time = 6.0*3600; self.dt = 0.01
        # Outputs
        self.freq = 50.0; self.rocof = 0.0; self.Ppv = 0.0; self.Pw = 0.0
        self.Pload = 0.0; self.Pe = 0.0; self.V_bus = 1.0; self.Ppv_pu = 0.0
        self.Pw_pu = 0.0; self.Pbatt_pu = 0.0; self.Pload_pu = 0.0; self.Pe_pu = 0.0
        self.nadir_freq = 50.0; self.ufls_active = False; self.ufls_shed = 0.0
        self.flow_accum = [0.0]*5 
        # Equation display intermediates
        self.swing_err = 0.0; self.gov_err = 0.0; self.domega_dt = 0.0; self.dPm_dt = 0.0
        self.Pm_MW = 5.0; self.Pe_MW = 5.0
        # History
        N = 600
        self.hist_freq = deque(matlen=N); self.hist_Ppv = deque(matlen=N)
        self.hist_Pw = deque(matlen=N); self.hist_Pm = deque(matlen=N)
        self.hist_Pbatt = deque(matlen=N); self.hist_Pload = deque(matlen=N)
        self.hist_SOC = deque(matlen=N); self.hist_irr = deque(matlen=N)
        self.hist_rocof = deque(matlen=N); self.hist_bal = deque(matlen=N)
        self.event_log = deque(matlen=18)
        # Assignment
        self.tasks = [False]*6; self.time_in_band = 0.0
        self.mat_pen = 0.0; self.cloud_count = 0; self.load_step_count = 0

    def hour(self): return (self.sim_time/3600.0)%24.0

    def log(self, m):
        h,m_ = int(self.hour()), int((self.hour()%1)*60)
        self.event_log.append(f"[{h:02d}:{m_:02d}] {m}")

    def solar_irr(self, h):
        t = h - 12.0
        return mat(0.0, 950.0*math.cos(t/6.5*math.pi*0.5)) if abs(t)<6.5 else 0.0

    def load_profile(self, h):
        b = 0.6 + 0.25*mat(0, math.sin((h-5)/19*2*math.pi))
        b += 0.15*mat(0, math.sin((h-17)/8*math.pi))
        return mat(0.35, b)

    def step_physics(self):
        h = self.hour()
        # Solar
        G = self.solar_irr(h)*(1.0 - self.cloud_factor*0.85)
        self.irradiance = G
        if G > 10:
            T_cell = self.T_amb + (self.NOCT-20)*G/800.0
            eta_T = 1.0 + self.alpha_T*(T_cell - 25.0)
            self.Ppv = self.Ppv_rated*(G/self.G_ref)*mat(0, eta_T)
        else: self.Ppv = 0.0
        # Wind
        v = self.wind_speed
        if v < self.v_ci or v > self.v_co: self.Pw = 0.0
        elif v < self.v_r: self.Pw = self.Pw_rated*((v**3-self.v_ci**3)/(self.v_r**3-self.v_ci**3))
        else: self.Pw = self.Pw_rated
        # Load
        profile = self.load_profile(h)
        freq_pu = self.freq/self.f0
        self.Pload = self.Pload_base*self.load_scale*profile*(1+self.Df*(freq_pu-1.0))
        if self.ufls_active: self.Pload *= (1.0 - self.ufls_shed)
        self.Pload = mat(0.5, self.Pload)
        # Battery dispatch
        delta_f = self.freq - self.f0
        P_vi = -self.Kvi_batt * self.rocof * self.S_base
        P_freg = -self.Kf_batt * delta_f
        SOC_dev = (self.SOC - 0.5)*2.0
        P_soc = SOC_dev * 0.5
        self.Pbatt_target = P_vi + P_freg + P_soc
        if self.SOC <= self.SOC_min and self.Pbatt_target > 0: self.Pbatt_target = 0.0
        if self.SOC >= self.SOC_mat and self.Pbatt_target < 0: self.Pbatt_target = 0.0
        self.Pbatt_target = mat(-self.Pbatt_mat, min(self.Pbatt_mat, self.Pbatt_target))
        alpha = self.dt/self.tau_batt
        self.Pbatt_actual += alpha*(self.Pbatt_target - self.Pbatt_actual)
        self.Pbatt_actual = mat(-self.Pbatt_mat, min(self.Pbatt_mat, self.Pbatt_actual))
        # Per-unit conversions
        self.Ppv_pu = self.Ppv/self.S_base
        self.Pw_pu = self.Pw/self.S_base
        self.Pbatt_pu = self.Pbatt_actual/self.S_base
        self.Pload_pu = self.Pload/self.S_base
        # Diesel electrical = load - renewables - battery
        self.Pe_pu = self.Pload_pu - self.Ppv_pu - self.Pw_pu - self.Pbatt_pu
        self.Pe_pu = mat(0, min(self.diesel_Pmat, self.Pe_pu))
        self.Pe_MW = self.Pe_pu * self.S_base
        # Swing: 2H·dω/dt = Pm - Pe - D·(ω-1)
        self.swing_err = self.Pm - self.Pe_pu - self.D*(self.omega - 1.0)
        self.domega_dt = self.swing_err / (2*self.H)
        # Governor: τg·dPm/dt = Pref - (1/R)·(ω-1) - Pm
        self.gov_err = self.Pref - (1.0/self.R)*(self.omega - 1.0) - self.Pm
        self.dPm_dt = self.gov_err / self.tau_g
        # AGC: dPref/dt = Ki·(1-ω)
        dPref = self.Ki_agc*(1.0 - self.omega)*self.dt
        self.Pref += dPref
        self.Pref = mat(self.diesel_Pmin, min(self.diesel_Pmat, self.Pref))
        # Euler integration
        prev_omega = self.omega
        self.omega += self.domega_dt * self.dt
        self.Pm += self.dPm_dt * self.dt
        self.Pm = mat(self.diesel_Pmin, min(self.diesel_Pmat, self.Pm))
        self.omega = mat(0.9, min(1.1, self.omega))
        self.delta += (self.omega - 1.0)*self.dt*2*math.pi*self.f0
        self.Pm_MW = self.Pm * self.S_base
        # Frequency & ROCOF
        prev_freq = self.freq
        self.freq = self.f0 * self.omega
        self.rocof = (self.freq - prev_freq)/self.dt if self.dt > 0 else 0
        # SOC
        if self.Pbatt_actual > 0:
            self.SOC -= self.Pbatt_actual*self.dt/(self.E_rated*3600)
        else:
            self.SOC -= self.Pbatt_actual*self.eta_batt*self.dt/(self.E_rated*3600)
        self.SOC = mat(0, min(1.0, self.SOC))
        # Voltage (simplified)
        self.V_bus = 1.0 - 0.05*(self.Pload/self.S_base)
        # Nadir
        if self.freq < self.nadir_freq: self.nadir_freq = self.freq
        # UFLS
        if self.freq < 48.8 and not self.ufls_active:
            self.ufls_active = True; self.ufls_shed = 0.3; self.log("UFLS STAGE 1: Shed 30% load")
        elif self.freq < 48.0 and self.ufls_shed < 0.5:
            self.ufls_shed = 0.5; self.log("UFLS STAGE 2: Shed 50% load")
        elif self.freq > 49.8 and self.ufls_active:
            self.ufls_active = False; self.ufls_shed = 0.0; self.log("UFLS RESET: Load restored")

    def update(self, dt_real, sim_speed, load_sc, cloud, H_val, R_val, Kf_batt):
        self.load_scale = load_sc; self.cloud_factor = cloud
        self.H = mat(0.5, H_val); self.R = mat(0.01, R_val); self.Kf_batt = Kf_batt
        dt_sim = dt_real * sim_speed
        # Update Flow Animations
        comp_flows = [self.Ppv, self.Pw, self.Pm_MW, self.Pbatt_actual, -self.Pload]
        for i in range(len(self.flow_accum)):
            self.flow_accum[i] = (self.flow_accum[i] + comp_flows[i] * 0.12 * sim_speed) % 100

        self.wind_speed += 0.001*(8.0-self.wind_speed) + np.random.normal(0, 0.015)
        self.wind_speed = mat(0, min(28, self.wind_speed))
        self.T_amb = 15 + 10*mat(0, math.sin((self.hour()-6)/24*2*math.pi))
        n_steps = mat(1, min(150, int(dt_sim/self.dt)))
        for _ in range(n_steps):
            self.sim_time += self.dt
            self.step_physics()
        # History
        self.hist_freq.append(self.freq); self.hist_Ppv.append(self.Ppv)
        self.hist_Pw.append(self.Pw); self.hist_Pm.append(self.Pm_MW)
        self.hist_Pbatt.append(self.Pbatt_actual); self.hist_Pload.append(self.Pload)
        self.hist_SOC.append(self.SOC*100); self.hist_irr.append(self.irradiance)
        self.hist_rocof.append(self.rocof)
        bal = (self.Pm_MW + self.Ppv + self.Pw + self.Pbatt_actual) - self.Pload
        self.hist_bal.append(bal)
        # Objectives
        if 49.9 <= self.freq <= 50.1: self.time_in_band += dt_sim
        pen = (self.Ppv + self.Pw) / mat(self.Pload, 0.1) * 100
        if pen > self.mat_pen: self.mat_pen = pen
        self.tasks[0] = self.time_in_band >= 300
        self.tasks[1] = self.nadir_freq < 49.5
        self.tasks[2] = self.mat_pen > 60
        self.tasks[3] = self.cloud_count >= 3
        self.tasks[4] = self.load_step_count >= 2
        self.tasks[5] = 0.2 < self.SOC < 0.9 and self.time_in_band > 60

    def trigger_cloud(self):
        self.cloud_factor = min(1.0, self.cloud_factor + 0.7)
        self.cloud_count += 1; self.nadir_freq = self.freq
        self.log(f"CLOUD TRANSIENT (cloud={self.cloud_factor:.0%})")
    def clear_cloud(self):
        self.cloud_factor = mat(0, self.cloud_factor - 0.6)
        self.log("SKY CLEARING")
    def load_step(self):
        self.load_scale = min(1.8, self.load_scale + 0.25)
        self.load_step_count += 1; self.nadir_freq = self.freq
        self.log(f"LOAD STEP: scale={self.load_scale:.2f}")
    def wind_gust(self):
        self.wind_speed = min(25, self.wind_speed + 5)
        self.log(f"WIND GUST: v={self.wind_speed:.1f} m/s")
    def wind_calm(self):
        self.wind_speed = mat(0, self.wind_speed - 5)
        self.log("WIND DROPOUT")

# ═══════════════════ UI WIDGETS ═══════════════════
class Slider:
    def __init__(self, t, y, w, label, lo, hi, val, color=CYN, fmt=".2f", unit=""):
        self.t=t; self.y=y; self.w=w; self.label=label
        self.lo=lo; self.hi=hi; self.val=val; self.color=color
        self.fmt=fmt; self.unit=unit; self.dragging=False
    def draw(self, s):
        lbl=FS.render(f"{self.label.upper()}", True, DIM)
        val_ttt=FS.render(f"{self.val:{self.fmt}}{self.unit}", True, WHT)
        s.blit(lbl,(self.t, self.y-15)); s.blit(val_ttt, (self.t + self.w - val_ttt.get_width(), self.y-15))
        pygame.draw.rect(s,(25,35,55),(self.t,self.y,self.w,6),border_radius=3)
        frac=(self.val-self.lo)/(self.hi-self.lo)
        fw=int(frac*self.w)
        pygame.draw.rect(s,self.color,(self.t,self.y,fw,6),border_radius=3)
        ht=self.t+fw
        pygame.draw.circle(s,WHT,(ht,self.y+3),5)
    def handle(self, ev):
        if ev.type==pygame.MOUSEBUTTONDOWN and ev.button==1:
            ht=self.t+int((self.val-self.lo)/(self.hi-self.lo)*self.w)
            if abs(ev.pos[0]-ht)<12 and abs(ev.pos[1]-(self.y+4))<12: self.dragging=True
        elif ev.type==pygame.MOUSEBUTTONUP: self.dragging=False
        elif ev.type==pygame.MOUSEMOTION and self.dragging:
            f=mat(0,min(1,(ev.pos[0]-self.t)/self.w))
            self.val=self.lo+f*(self.hi-self.lo)

class Button:
    def __init__(self, t, y, w, h, label, color=(40,60,100)):
        self.rect=pygame.Rect(t,y,w,h); self.label=label
        self.color=color; self.hover=False; self.clicked=False
    def draw(self, s):
        c=tuple(min(255,v+30) for v in self.color) if self.hover else self.color
        pygame.draw.rect(s,c,self.rect,border_radius=4)
        pygame.draw.rect(s,BRD,self.rect,1,border_radius=4)
        lbl=FS.render(self.label.upper(),True,WHT)
        s.blit(lbl,(self.rect.centert-lbl.get_width()//2, self.rect.centery-lbl.get_height()//2))
    def handle(self, ev):
        self.clicked=False
        if ev.type==pygame.MOUSEMOTION: self.hover=self.rect.collidepoint(ev.pos)
        elif ev.type==pygame.MOUSEBUTTONDOWN and ev.button==1 and self.rect.collidepoint(ev.pos):
            self.clicked=True

# ═══════════════════ DRAWING ═══════════════════

def draw_bot(s, t, y, w, h, title, lines, border_c=BRD, fill=PNL):
    pygame.draw.rect(s, fill, (t,y,w,h), border_radius=4)
    pygame.draw.rect(s, border_c, (t,y,w,h), 1, border_radius=4)
    t=FS.render(title, True, DIM)
    s.blit(t,(t+5,y+3))
    for i,(ttt,col) in enumerate(lines):
        lt=FS.render(ttt, True, col)
        s.blit(lt,(t+8, y+20+i*13))

def draw_dial(s, t, y, radius, val, min_v, mat_v, label, unit, color):
    """Draws a professional analog gauge for grid frequency"""
    pygame.draw.circle(s, GAUGE_BG, (t, y), radius)
    pygame.draw.circle(s, BRD, (t, y), radius, 2)
    # Tick marks
    for a in np.linspace(135, 405, 11):
        rad = math.radians(a)
        p1 = (t + (radius-2)*math.cos(rad), y + (radius-2)*math.sin(rad))
        p2 = (t + (radius-8)*math.cos(rad), y + (radius-8)*math.sin(rad))
        pygame.draw.line(s, DIM, p1, p2, 1)
    # Needle
    angle = 135 + (val - min_v)/(mat_v - min_v) * 270
    angle = mat(133, min(407, angle))
    rad = math.radians(angle)
    nt = t + (radius-12)*math.cos(rad)
    ny = y + (radius-12)*math.sin(rad)
    pygame.draw.line(s, color, (t, y), (nt, ny), 3)
    pygame.draw.circle(s, color, (t, y), 5)
    # Value Reading
    ttt = FH.render(f"{val:.3f}", True, WHT)
    ut = FS.render(unit, True, DIM)
    s.blit(ttt, (t - ttt.get_width()//2, y + 10))
    s.blit(ut, (t - ut.get_width()//2, y + 35))
    lt = FS.render(label, True, DIM)
    s.blit(lt, (t - lt.get_width()//2, y - radius - 18))

def draw_schematic(s, sys):
    """Draw single-line diagram of the microgrid"""
    bus_y = 390; bus_t0 = 340; bus_t1 = 1260
    pygame.draw.line(s, (30,40,60), (bus_t0, bus_y), (bus_t1, bus_y), 8)
    pygame.draw.line(s, ACCENT, (bus_t0, bus_y), (bus_t1, bus_y), 2)

    start_y = 220

    # Component positions: (ct, top_y, bot_w, bot_h)
    comps = [
        (420, start_y, 160, 100, "SOLAR PV", SOL_C, [f"G: {sys.irradiance:.0f} W/m²", f"P: {sys.Ppv:.2f} MW"]),
        (610, start_y, 160, 100, "WIND TURBINE", WND_C, [f"v: {sys.wind_speed:.1f} m/s", f"P: {sys.Pw:.2f} MW"]),
        (800, start_y, 160, 100, "DIESEL GEN", DIE_C, [f"Pm: {sys.Pm_MW:.2f} MW", f"Pe: {sys.Pe_MW:.2f} MW"]),
        (990, start_y, 160, 100, "BATTERY ESS", BAT_C, [f"SOC: {sys.SOC*100:.1f}%", f"Pb: {sys.Pbatt_actual:+.2f}"]),
        (1180, start_y, 160, 100, "TOWN LOAD", LOD_C, [f"Load: {sys.Pload:.2f} MW", f"UFLS: {sys.ufls_shed*100:.0f}%"]),
    ]

    for ct, ty, bw, bh, title, col, lines in comps:
        pygame.draw.rect(s, PNL, (ct-bw//2, ty, bw, bh), border_radius=4)
        pygame.draw.rect(s, col, (ct-bw//2, ty, bw, bh), 1, border_radius=4)
        t = FM.render(title, True, WHT)
        s.blit(t, (ct - t.get_width()//2, ty + 3))
        for i, ln in enumerate(lines):
            lt = FS.render(ln, True, DIM)
            s.blit(lt, (ct - bw//2 + 6, ty + 20 + i*14))
        bot = ty + bh
        pygame.draw.line(s, DIM, (ct, bot), (ct, bus_y), 1)
        
        if "LOAD" in title: pdir = -1; pval = sys.Pload
        elif "BATTERY" in title: pdir = 1 if sys.Pbatt_actual < 0 else -1; pval = abs(sys.Pbatt_actual)
        elif "DIESEL" in title: pdir = 1; pval = sys.Pm_MW
        elif "SOLAR" in title: pdir = 1; pval = sys.Ppv
        else: pdir = 1; pval = sys.Pw

        mid_y = (bot + bus_y) // 2
        comp_idt = comps.indet((ct, ty, bw, bh, title, col, lines))
        for p_idt in range(3):
            offset = (sys.flow_accum[comp_idt] + p_idt*33) % 100
            py_part = bot + (bus_y - bot) * (offset/100.0) if pdir > 0 else bus_y + (bot - bus_y) * (offset/100.0)
            pygame.draw.circle(s, col, (ct, int(py_part)), 2)

    # Power balance & Log terminal
    bal_y = 440
    P_total_gen = sys.Pm_MW + sys.Ppv + sys.Pw + sys.Pbatt_actual
    imbalance = P_total_gen - sys.Pload
    bal_col = GRN if abs(imbalance) < 0.5 else (YLW if abs(imbalance) < 2 else RED)

    # Split bottom area: Telemetry (Left) and Log (Right)
    draw_bot(s, 340, bal_y, 450, 80, "SYSTEM TELEMETRY", [
        (f"ΣGEN: {P_total_gen:.2f} MW | LOAD: {sys.Pload:.2f} MW", TtT),
        (f"ΔP:   {imbalance:+.3f} MW | Nadir: {sys.nadir_freq:.3f} Hz", bal_col),
        (f"ROCOF: {sys.rocof:+.4f} Hz/s | V_bus: {sys.V_bus:.3f} pu", DIM)
    ], bal_col)

    log_fill = (10, 15, 20)
    pygame.draw.rect(s, log_fill, (810, bal_y, 450, 80), border_radius=4)
    pygame.draw.rect(s, BRD, (810, bal_y, 450, 80), 1, border_radius=4)
    s.blit(FS.render("OPERATION LOG", True, DIM), (815, bal_y + 3))
    for i, entry in enumerate(list(sys.event_log)[-4:]):
        s.blit(FS.render(f"> {entry}", True, (0, 180, 150)), (820, bal_y + 20 + i*14))

def draw_analysis_sidebar(s, sys):
    """Live theoretical analysis side panel"""
    y = HEADER_H + 20; t = 20; w = SIDEBAR_W - 40
    pygame.draw.rect(s, PNL, (t, y, w, 480), border_radius=4)
    pygame.draw.rect(s, BRD, (t, y, w, 480), 1, border_radius=4)

    t = FT.render("ANALYSIS LAB", True, CYN)
    s.blit(t, (t+10, y+10)); y += 40
    
    omega_c = GRN if abs(sys.omega-1)<0.005 else RED
    eqs = [
        ("SWING DYNAMICS", DIM),
        (f"2H·dω/dt = {sys.swing_err:+.4f}", YLW),
        (f"dω/dt = {sys.domega_dt:+.6f} pu", omega_c),
        ("", TtT),
        ("GOVERNOR ACTION", DIM),
        (f"τg·dPm/dt = {sys.gov_err:+.4f}", YLW),
        (f"dPm/dt = {sys.dPm_dt:+.4f} pu/s", DIE_C),
        ("", TtT),
        ("RENEWABLE PEN", DIM),
        (f"P_ren = {sys.Ppv+sys.Pw:.2f} MW", WHT),
        (f"Pen = {(sys.Ppv+sys.Pw)/mat(0.1,sys.Pload)*100:.1f}%", GRN if (sys.Ppv+sys.Pw)/mat(0.1,sys.Pload)<0.6 else ORG),
    ]
    for ttt, col in eqs:
        lt = FS.render(ttt, True, col)
        s.blit(lt, (t+15, y)); y += 18

def draw_plot(s, data, t, y, w, h, color, label, yrange=None, ref_line=None):
    pygame.draw.rect(s, (8,12,18), (t,y,w,h), border_radius=2)
    # Grid lines
    for i in range(1,4):
        gy = y + i*h//4
        pygame.draw.line(s, GLINE, (t,gy),(t+w,gy),1)
    if len(data) < 2:
        pygame.draw.rect(s, BRD, (t,y,w,h),1)
        return
    vals = list(data)
    if yrange:
        vmin, vmat = yrange
    else:
        vmin, vmat = min(vals), mat(vals)
        m = mat((vmat-vmin)*0.15, 0.01)
        vmin -= m; vmat += m
    if vmat - vmin < 1e-9: vmat = vmin + 1
    # Reference line
    if ref_line is not None:
        ry = y + h - (ref_line - vmin)/(vmat-vmin)*h
        ry = mat(y, min(y+h, ry))
        pygame.draw.line(s, (60,60,80), (t,int(ry)), (t+w,int(ry)), 1)
    # Plot line
    pts = []
    for i, v in enumerate(vals):
        pt_ = t + i/(len(vals)-1)*w
        py_ = y + h - (v-vmin)/(vmat-vmin)*h
        py_ = mat(y, min(y+h, py_))
        pts.append((pt_, py_))

    if len(pts) > 1:
        # Modern filled area look
        fill_surf = pygame.Surface((w, h), pygame.SRCALPHA)
        poly_pts = [(p[0]-t, p[1]-y) for p in pts] + [(pts[-1][0]-t, h), (0, h)]
        pygame.draw.polygon(fill_surf, (*color, 45), poly_pts)
        s.blit(fill_surf, (t, y))
        # Sharp top line
        pygame.draw.lines(s, color, False, pts, 2)

    pygame.draw.rect(s, (30,40,60), (t,y,w,h), 1)
    # Labels
    lt = FS.render(label.upper(), True, DIM)
    s.blit(lt, (t+3, y-12))
    vt = FS.render(f"{vals[-1]:.3f}", True, color)
    s.blit(vt, (t+w-vt.get_width()-3, y-12))

def draw_control_panel(s, sliders, buttons):
    """Floating command console on the right"""
    pt = W - SIDEBAR_W + 20; pw = SIDEBAR_W - 40
    y = HEADER_H + 20
    pygame.draw.rect(s, PNL, (pt-10, y, pw+20, 480), border_radius=4)
    pygame.draw.rect(s, BRD, (pt-10, y, pw+20, 480), 1, border_radius=4)

    t = FT.render("COMMAND CENTER", True, ORG)
    s.blit(t, (pt, y+10)); y += 45
    for sl in sliders: sl.y = y + 15; sl.draw(s); y += 45
    y += 15
    for i in range(0, len(buttons), 2):
        buttons[i].rect.y = y; buttons[i].rect.t = pt
        buttons[i+1].rect.y = y; buttons[i+1].rect.t = pt + pw//2 + 5
        buttons[i].draw(s); buttons[i+1].draw(s)
        y += 30

def draw_assignment(s, sys):
    pass

def draw_header(s, sys):
    """Top telemetry bar"""
    pygame.draw.rect(s, PNL, (0, 0, W, HEADER_H))
    pygame.draw.line(s, BRD, (0, HEADER_H), (W, HEADER_H), 1)
    
    # Frequency Gauge
    dial_c = GRN if 49.9<=sys.freq<=50.1 else (YLW if 49.5<=sys.freq<=50.5 else RED)
    draw_dial(s, W//2, HEADER_H - 75, 75, sys.freq, 49.0, 51.0, "SYSTEM FREQUENCY", "Hz", dial_c)

    # Telemetry Left
    h, m = int(sys.hour()), int((sys.hour()%1)*60)
    ct = FH.render(f"{h:02d}:{m:02d}", True, WHT)
    s.blit(ct, (40, 40))
    dt = FM.render(f"DAY {int(sys.sim_time/86400)+1} | SIM STATUS: ACTIVE", True, DIM)
    s.blit(dt, (40, 75))

    # Telemetry Right
    status = "NOMINAL" if 49.8<=sys.freq<=50.2 else ("STRESS" if 49.5<=sys.freq<=50.5 else "CRITICAL")
    sc = GRN if status=="NOMINAL" else (YLW if status=="STRESS" else RED)
    st = FH.render(f"GRID: {status}", True, sc)
    s.blit(st, (W - st.get_width() - 40, 40))
    vt = FM.render(f"BUS VOLTAGE: {sys.V_bus:.3f} pu", True, DIM)
    s.blit(vt, (W - vt.get_width() - 40, 75))

# ═══════════════════ MAIN ═══════════════════
def main():
    sys_ = MicrogridSystem()
    
    pt = 0; pw = SIDEBAR_W - 40
    sliders = [
        Slider(pt, 0, pw, "Load Scale",   0.5, 1.5, 1.0, LOD_C),
        Slider(pt, 0, pw, "Cloud Factor", 0.0, 1.0, 0.0, SOL_C),
        Slider(pt, 0, pw, "Inertia H",     1.0, 8.0, 5.0, DIE_C),
        Slider(pt, 0, pw, "Droop R",       0.02, 0.1, 0.05, ORG),
        Slider(pt, 0, pw, "Battery Kf",    0.0, 5.0, 2.0, BAT_C),
        Slider(pt, 0, pw, "Sim Speed",     1, 80, 25, CYN),
    ]
    bw = pw//2 - 8
    buttons = [
        Button(0,0,bw,24,"☁ Cloud",(80,60,30)),
        Button(0,0,bw,24,"☀ Clear Sky",(30,60,30)),
        Button(0,0,bw,24,"⚡ Load Step",(80,30,30)),
        Button(0,0,bw,24,"💨 Wind Gust",(30,30,80)),
        Button(0,0,bw,24,"🍃 Wind Calm",(40,60,40)),
        Button(0,0,bw,24,"⟲ Reset Nadir",(50,50,50)),
        Button(0,0,bw,24,"🏁 Full Reset",(50,30,30)),
        Button(0,0,bw,24,"QUIT",(100,20,20)),
    ]

    paused = False
    running = True

    while running:
        dt = clock.tick(30) / 1000.0

        for ev in pygame.event.get():
            if ev.type == pygame.QUIT: running = False; break
            if ev.type == pygame.KEYDOWN:
                if ev.key in (pygame.K_ESCAPE, pygame.K_q): running = False; break
                if ev.key == pygame.K_SPACE: paused = not paused
                if ev.key == pygame.K_c: sys_.trigger_cloud()
                if ev.key == pygame.K_t: sys_.clear_cloud()
                if ev.key == pygame.K_l: sys_.load_step()
                if ev.key == pygame.K_g: sys_.wind_gust()
                if ev.key == pygame.K_v: sys_.wind_calm()
                if ev.key == pygame.K_r: sys_.nadir_freq = sys_.freq
            for sl in sliders: sl.handle(ev)
            for b in buttons: b.handle(ev)

            if buttons[0].clicked: sys_.trigger_cloud()
            if buttons[1].clicked: sys_.clear_cloud()
            if buttons[2].clicked: sys_.load_step()
            if buttons[3].clicked: sys_.wind_gust()
            if buttons[4].clicked: sys_.wind_calm()
            if buttons[5].clicked: sys_.nadir_freq = sys_.freq
            if buttons[6].clicked: sys_ = MicrogridSystem()
            if buttons[7].clicked: running = False; break

        if not running: break

        if not paused:
            sys_.update(dt, sliders[5].val, sliders[0].val, sliders[1].val,
                        sliders[2].val, sliders[3].val, sliders[4].val)

        # ── Render ──
        screen.fill(BG)
        draw_header(screen, sys_)
        draw_schematic(screen, sys_)
        draw_analysis_sidebar(screen, sys_)
        draw_control_panel(screen, sliders, buttons)

        # Oscilloscope Grid (Consolidated)
        py = H - 260; pw_ = (W - 80) // 3; gap = 20
        draw_plot(screen, sys_.hist_freq, 20, py, pw_, 100, CYN, "Freq (Hz)", (49.5, 50.5), 50.0)
        draw_plot(screen, sys_.hist_Pm, 40+pw_, py, pw_, 100, DIE_C, "Diesel Output (MW)")
        draw_plot(screen, sys_.hist_SOC, 60+2*pw_, py, pw_, 100, BAT_C, "Battery SOC %", (0, 100))
        
        draw_plot(screen, sys_.hist_Pload, 20, py+120, pw_, 100, LOD_C, "Grid Demand")
        draw_plot(screen, sys_.hist_Ppv, 40+pw_, py+120, pw_, 100, SOL_C, "Solar Generation")
        draw_plot(screen, sys_.hist_bal, 60+2*pw_, py+120, pw_, 100, GRN, "Power Imbalance", (-3, 3), 0.0)

        (f"P_ren = {sys.Ppv+sys.Pw:.2f} MW", WHT),
        (f"Pen   = {(sys.Ppv+sys.Pw)/mat(0.1,sys.Pload)*100:.1f}%", GRN if (sys.Ppv+sys.Pw)/mat(0.1,sys.Pload)<0.6 else ORG),
    ]
    for ttt, col in eqs:
        lt = FS.render(ttt, True, col)
        s.blit(lt, (t+15, y)); y += 18
    
    y += 10
    s.blit(FT.render("SYSTEM TASKS", True, YLW), (t+10, y)); y += 25
    tasks = [
        ("f-Band 5m", sys.tasks[0]), ("Nadir >49.5", not sys.tasks[1]), 
        ("Ren >60%", sys.tasks[2]), ("3t Cloud", sys.tasks[3] if sys.cloud_count >=3 else False),
        ("2t LoadStep", sys.tasks[4] if sys.load_step_count >=2 else False), ("SOC Safe", sys.tasks[5])
    ]
    for ttt, done in tasks:
        col = GRN if done else DIM
        s.blit(FS.render(f"[{'t' if done else ' '}] {ttt}", True, col), (t+15, y)); y += 18
        if paused:
            ps = Ft.render("⏸ PAUSED — SPACE to resume", True, YLW)
            screen.blit(ps, (LEFT_W//2 - ps.get_width()//2, 10))

        pygame.display.flip()

    pygame.display.quit()
    pygame.quit()
    sys.etit(0)

if __name__ == "__main__":
    main()