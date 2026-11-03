"""
02_shear_and_density.py
- Compute wind shear exponent (power law) from concurrent multi-height data.
- Compute actual air density from station pressure/temperature (elevation 3040 m).
- Extrapolate 40 m wind speed series to hub height using derived shear exponent.
"""
import pandas as pd
import numpy as np

df = pd.read_parquet(r"D:\Aaryan\College\College Subjects\Wind Energy Conversion System\Lab\Tangbe wind data.parquet")

heights = np.array([20, 40, 60, 80])
speed_cols = [f"a{h}_wind_speed_mean" for h in heights]

# --- Wind shear exponent (power law): v2/v1 = (h2/h1)^alpha ---
# Use timestamps where ALL FOUR heights are valid simultaneously (concurrent fit),
# fit ln(v) = ln(v_ref) + alpha * ln(h/h_ref) via least squares per timestamp.
sub = df[speed_cols].dropna()
print(f"Timestamps with all 4 heights valid concurrently: {len(sub)} "
      f"({100*len(sub)/len(df):.1f}% of records)")

log_h = np.log(heights)
alphas = []
for _, row in sub.iterrows():
    v = row.values.astype(float)
    v = np.where(v <= 0, np.nan, v)
    if np.isnan(v).any():
        continue
    log_v = np.log(v)
    # linear regression log_v = a + alpha*log_h
    A = np.vstack([np.ones_like(log_h), log_h]).T
    coef, *_ = np.linalg.lstsq(A, log_v, rcond=None)
    alphas.append(coef[1])

alphas = np.array(alphas)
alpha_mean = np.nanmean(alphas)
alpha_median = np.nanmedian(alphas)
print(f"\nShear exponent alpha: mean={alpha_mean:.4f}, median={alpha_median:.4f}, "
      f"std={np.nanstd(alphas):.4f}  (n={len(alphas)})")

# Sanity check using simple 40m-80m pair over their common valid period
pair = df[["a40_wind_speed_mean", "a80_wind_speed_mean"]].dropna()
alpha_pair = np.log(pair["a80_wind_speed_mean"] / pair["a40_wind_speed_mean"]) / np.log(80/40)
alpha_pair = alpha_pair.replace([np.inf, -np.inf], np.nan).dropna()
print(f"Cross-check (40m-80m pair only): alpha median = {alpha_pair.median():.4f}  (n={len(alpha_pair)})")

ALPHA = round(float(alpha_median), 3)
print(f"\n--> Using alpha = {ALPHA} for extrapolation")

# --- Air density from actual station data (ideal gas law) ---
# rho = p / (R_specific * T), p in Pa, T in Kelvin, R_specific(dry air) = 287.05 J/(kg K)
R = 287.05
p_hpa = df["p5_air_pressure_mean"]          # hPa
t_c = df["t5_temperature_mean"]              # deg C
rho = (p_hpa * 100) / (R * (t_c + 273.15))   # kg/m^3

print(f"\n=== Air density (from station pressure/temp @ 5 m, elevation 3040 m) ===")
print(f"  mean rho = {rho.mean():.4f} kg/m3")
print(f"  min/max  = {rho.min():.4f} / {rho.max():.4f} kg/m3")
print(f"  (sea-level standard = 1.225 kg/m3 -> site density is "
      f"{100*rho.mean()/1.225:.1f}% of sea level)")

RHO_MEAN = float(rho.mean())

# --- Extrapolate 40 m series (most complete) to hub height using ALPHA ---
HUB_HEIGHT = 100  # m -- chosen hub height for candidate turbine (see script 03 for rationale)
v40 = df["a40_wind_speed_mean"]
v_hub = v40 * (HUB_HEIGHT / 40) ** ALPHA

print(f"\n=== Extrapolated hub-height ({HUB_HEIGHT} m) wind speed from 40 m series ===")
print(f"  n valid = {v_hub.notna().sum()}")
print(f"  mean = {v_hub.mean():.3f} m/s, std = {v_hub.std():.3f} m/s")
print(f"  min/max = {v_hub.min():.3f} / {v_hub.max():.3f} m/s")

out = pd.DataFrame({
    "v40": v40,
    "v_hub": v_hub,
    "rho": rho,
})
out.attrs["alpha"] = ALPHA
out.attrs["rho_mean"] = RHO_MEAN
out.attrs["hub_height"] = HUB_HEIGHT
out.to_parquet("/home/claude/wecs_project/data/hub_series.parquet")

# save scalar summary for later scripts
import json
summary = {"alpha": ALPHA, "rho_mean": RHO_MEAN, "hub_height": HUB_HEIGHT,
           "elevation_m": 3040, "n_hub_valid": int(v_hub.notna().sum())}
with open("/home/claude/wecs_project/data/site_constants.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\nSaved hub_series.parquet and site_constants.json")