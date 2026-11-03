import pandas as pd
from windrose import WindroseAxes
import matplotlib.pyplot as plt

# Use the exact name of the sheet instead of the number 3
df = pd.read_excel(r'D:\Aaryan\College\College Subjects\Wind Energy Conversion System\Lab\Tangbe wind data.xlsx', sheet_name='10Tangbe')

# Drop rows where wind direction or wind speed is missing (NaN)
df_clean = df.dropna(subset=['d78.5_wind_direction_mean', 'a80_wind_speed_mean'])

# Create the Windrose plot
ax = WindroseAxes.from_ax()
ax.bar(df_clean['d78.5_wind_direction_mean'], 
       df_clean['a80_wind_speed_mean'], 
       normed=True, 
       opening=0.9, 
       edgecolor='white', 
       bins=[0, 2, 4, 6, 8, 10, 12])

ax.set_legend()
plt.tight_layout() # Prevents the legend from being cut off
plt.show()