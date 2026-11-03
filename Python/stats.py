# -*- coding: utf-8 -*-
"""
Expense Analysis & Forecasting for Excel File (by position)
-----------------------------------------------------------
Reads an Excel file with expenses from Feb 13, 2026 to Mar 6, 2026,
selects rows by position, cleans the data, and provides insights.

File path: D:\Aaryan\Drive\Personal\Expensese since feb 13 2026 or falgun 1 2082.xlsx
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# 1. CONFIGURATION – UPDATE PATH IF NECESSARY
# ============================================================
file_path = r"D:\Aaryan\Drive\Personal\Expensese since feb 13 2026 or falgun 1 2082.xlsx"

# ============================================================
# 2. LOAD DATA AND SELECT ROWS/COLUMNS BY POSITION
# ============================================================
# Read Excel, first column as index (Item names)
df_raw = pd.read_excel(file_path, index_col=0, header=0, engine='openpyxl')

# DEBUG: Print the first few index values to see what's there
print("First 15 index values from the file:")
print(df_raw.index[:15].tolist())
print("\n" + "="*60 + "\n")

# Select rows 1 through 10 (0‑based: rows 1 to 10 inclusive = iloc[1:11])
# and columns 0 through 19 (the first 20 date columns B:U)
df_cat = df_raw.iloc[1:11, :20].copy()

# Optional: reset index to clean category names (strip whitespace)
df_cat.index = df_cat.index.str.strip()

# Convert column names to dates (they may include time)
df_cat.columns = pd.to_datetime(df_cat.columns).date

print("Selected categories:")
print(df_cat.index.tolist())
print("\nNumber of date columns:", len(df_cat.columns))

# ============================================================
# 3. CLEAN EACH CELL: evaluate formulas and convert to numbers
# ============================================================
def clean_cell(val):
    """
    Convert a cell to a float.
    - If it's a string starting with '=', evaluate simple addition (e.g., "=50+20" → 70)
    - Otherwise try to convert to float directly.
    - Empty or non-convertible becomes 0.
    """
    if pd.isna(val):
        return 0.0
    if isinstance(val, str):
        val = val.strip()
        if val.startswith('='):
            # Evaluate only simple arithmetic (numbers and +)
            try:
                parts = val[1:].split('+')
                return sum(float(p.strip()) for p in parts)
            except:
                return 0.0
        else:
            try:
                return float(val)
            except:
                return 0.0
    else:
        return float(val)

# Apply cleaning to every cell
df_clean = df_cat.applymap(clean_cell)

# ============================================================
# 4. COMPUTE TOTALS
# ============================================================
daily_total = df_clean.sum(axis=0)          # sum over categories per day
category_total = df_clean.sum(axis=1)       # sum over days per category
total_expense = daily_total.sum()

print("=" * 60)
print("EXPENSE SUMMARY (Feb 13 – Mar 6, 2026)")
print("=" * 60)
print(f"Total expenses: Rs {total_expense:,.2f}")
print(f"Average daily expense: Rs {daily_total.mean():.2f}")
print(f"Max daily expense: Rs {daily_total.max():.2f} on {daily_total.idxmax()}")
print(f"Min daily expense: Rs {daily_total.min():.2f} on {daily_total.idxmin()}")
print("\nTop spending categories:")
print(category_total.sort_values(ascending=False).to_string(float_format="Rs %.2f"))

# ============================================================
# 5. VISUALIZATIONS
# ============================================================
plt.style.use('seaborn-v0_8-darkgrid')
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Daily total line chart
axes[0, 0].plot(daily_total.index, daily_total.values, marker='o', linestyle='-', color='b')
axes[0, 0].set_title('Daily Total Expenses')
axes[0, 0].set_xlabel('Date')
axes[0, 0].set_ylabel('Amount (Rs)')
axes[0, 0].tick_params(axis='x', rotation=45)

# Category pie chart
axes[0, 1].pie(category_total, labels=category_total.index, autopct='%1.1f%%', startangle=90)
axes[0, 1].set_title('Expense Breakdown by Category')

# Stacked area chart of categories over time
df_clean.T.plot.area(ax=axes[1, 0], legend=False)
axes[1, 0].set_title('Category Contributions Over Time')
axes[1, 0].set_xlabel('Date')
axes[1, 0].set_ylabel('Daily Amount (Rs)')
axes[1, 0].tick_params(axis='x', rotation=45)

# Daily total with 7‑day moving average
moving_avg = daily_total.rolling(window=7, min_periods=1).mean()
axes[1, 1].plot(daily_total.index, daily_total.values, 'o-', alpha=0.5, label='Daily total')
axes[1, 1].plot(daily_total.index, moving_avg.values, 'r-', linewidth=2, label='7‑day moving avg')
axes[1, 1].set_title('Daily Total with Moving Average')
axes[1, 1].set_xlabel('Date')
axes[1, 1].set_ylabel('Amount (Rs)')
axes[1, 1].legend()
axes[1, 1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.show()

# ============================================================
# 6. PREDICTION MODEL + 7-DAY FORECAST
# ============================================================
ts = daily_total.copy().astype(float)

def build_lag_features(series, max_lag=3):
    """
    Build supervised-learning features from a time series.
    Features:
    - lag_1, lag_2, ..., lag_n
    - 3-day rolling mean of previous days
    - day-of-week cyclical features (sin/cos)
    """
    df_feat = pd.DataFrame({'y': series})
    for lag in range(1, max_lag + 1):
        df_feat[f'lag_{lag}'] = df_feat['y'].shift(lag)
    df_feat['roll3'] = df_feat['y'].shift(1).rolling(3).mean()
    dow = pd.to_datetime(df_feat.index).dayofweek
    df_feat['dow_sin'] = np.sin(2 * np.pi * dow / 7)
    df_feat['dow_cos'] = np.cos(2 * np.pi * dow / 7)
    return df_feat.dropna()

max_lag = 3
df_model = build_lag_features(ts, max_lag=max_lag)
feature_cols = [f'lag_{i}' for i in range(1, max_lag + 1)] + ['roll3', 'dow_sin', 'dow_cos']

X_all = df_model[feature_cols]
y_all = df_model['y']

# Hold out last 5 observed days for quick validation
test_size = min(5, max(1, len(df_model) // 4))
X_train, X_test = X_all.iloc[:-test_size], X_all.iloc[-test_size:]
y_train, y_test = y_all.iloc[:-test_size], y_all.iloc[-test_size:]

model = LinearRegression()
model.fit(X_train, y_train)
val_pred = model.predict(X_test)
mae = mean_absolute_error(y_test, val_pred)

# Retrain on full data before forecasting
model.fit(X_all, y_all)
approx_pred_all = model.predict(X_all)

last_day = ts.index[-1]
future_dates = [last_day + timedelta(days=i) for i in range(1, 8)]

# Recursive multi-step forecast
history = ts.tolist()
future_pred_model = []
for fd in future_dates:
    lag_vals = [history[-i] for i in range(1, max_lag + 1)]
    roll3 = np.mean(history[-3:])
    dow = pd.to_datetime(fd).dayofweek
    row = lag_vals + [roll3, np.sin(2 * np.pi * dow / 7), np.cos(2 * np.pi * dow / 7)]
    pred = float(model.predict(np.array(row).reshape(1, -1))[0])
    pred = max(0.0, pred)
    future_pred_model.append(pred)
    history.append(pred)

# Baseline trend model for comparison
days_since_start = (ts.index - ts.index[0]).days
trend_model = LinearRegression()
trend_model.fit(days_since_start.values.reshape(-1, 1), ts.values)
future_days = [(d - ts.index[0]).days for d in future_dates]
future_pred_trend = trend_model.predict(np.array(future_days).reshape(-1, 1))

# Plot model approximation on historical data
plt.figure(figsize=(10, 5))
plt.plot(ts.index, ts.values, 'o-', color='steelblue', alpha=0.6, label='Actual historical')
plt.plot(df_model.index, approx_pred_all, 'r-', linewidth=2, label='Model approximation (fit)')
plt.plot(X_test.index, val_pred, 'ko', markersize=7, label='Validation predictions')
plt.title('Model Approximation vs Actual (Historical)')
plt.xlabel('Date')
plt.ylabel('Amount (Rs)')
plt.legend()
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# Plot future forecasts
plt.figure(figsize=(10, 5))
plt.plot(ts.index, ts.values, 'bo-', label='Historical')
plt.plot(future_dates, future_pred_model, 'ro--', label='Lag-Feature Model Forecast')
plt.plot(future_dates, future_pred_trend, 'go--', label='Trend Baseline Forecast')
plt.title('Daily Expense Forecast (Next 7 Days)')
plt.xlabel('Date')
plt.ylabel('Amount (Rs)')
plt.legend()
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

print("\n" + "=" * 60)
print("PREDICTION MODEL PERFORMANCE")
print("=" * 60)
print(f"Validation MAE (last {test_size} days): Rs {mae:.2f}")

print("\n" + "=" * 60)
print("FORECAST FOR NEXT 7 DAYS (Mar 7 - Mar 13, 2026)")
print("=" * 60)
forecast_df = pd.DataFrame({
    'Date': future_dates,
    'Lag-Feature Model (Rs)': future_pred_model,
    'Trend Baseline (Rs)': future_pred_trend
})
print(forecast_df.to_string(index=False, float_format="%.2f"))

# ============================================================
# 7. ADDITIONAL INSIGHTS
# ============================================================
# Average spending by weekday
daily_with_dow = daily_total.to_frame(name='amount')
daily_with_dow['day_of_week'] = daily_with_dow.index.dayofweek
weekly_avg = daily_with_dow.groupby('day_of_week')['amount'].mean()
days_map = {0: 'Mon', 1: 'Tue', 2: 'Wed', 3: 'Thu', 4: 'Fri', 5: 'Sat', 6: 'Sun'}
weekly_avg.index = weekly_avg.index.map(days_map)

print("\nAverage spending by weekday:")
print(weekly_avg.to_string(float_format="Rs %.2f"))

# Category variability (standard deviation)
variability = df_clean.std(axis=1).sort_values(ascending=False)
print("\nMost variable categories (standard deviation):")
print(variability.to_string(float_format="Rs %.2f"))
