import numpy as np
import pandas as pd

# Number of days and number of river cells
num_days = 730  # 2 years of data
num_cells = 200

# Generate time series (daily)
time_series = pd.date_range(start='2024-01-01', periods=num_days, freq='D')

# Generate synthetic tidal range and salinity gradient data
np.random.seed(42)
tidal_range = 2 + 0.5 * np.sin(np.linspace(0, 2 * np.pi, num_days)) + 0.1 * np.random.randn(num_days)
salinity_gradient = 30 + 5 * np.sin(np.linspace(0, 2 * np.pi, num_days))[:, np.newaxis] + 0.5 * np.random.randn(num_days, num_cells)

# Combine the data into a DataFrame
data = pd.DataFrame({
    'Date': time_series,
    'Tidal_Range': tidal_range
})
for i in range(num_cells):
    data[f'Salinity_Cell_{i+1}'] = salinity_gradient[:, i]

# Save to CSV
data.to_csv('calibration_data.csv', index=False)
print("Calibration data generated and saved to 'calibration_data.csv'")